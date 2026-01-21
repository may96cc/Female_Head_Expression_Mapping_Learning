import pandas as pd
import numpy as np
from scipy import stats
import re
import os
import sys
import argparse

"""
QTL Mapping with Permutation Testing (FAST + CORRECT)
======================================================
Output: F-stat, LOD, PVE (effect size)
Significance thresholds determined by permutation of LOD scores (Churchill & Doerge 1994)

Key insight: Vectorize permutations PER BLOCK using hat matrix trick,
then take MAX across blocks for each permutation to get genome-wide threshold.

Formulas (from Broman & Sen 2009, King & Long 2017):
- LOD from F-stat: LOD = (n/2) * log10(1 + F * df_reg / df_res)
- PVE from LOD: PVE = 1 - 10^(-2*LOD/n)

Permutation threshold approach (Churchill & Doerge 1994):
- For each permutation: permute phenotype, scan genome, record MAX LOD
- 95th percentile of max LOD distribution = 5% FWER threshold
"""


def f_to_lod(f_stat, n, df_reg, df_res):
    """
    Convert F-statistic to LOD score.
    
    Formula from Broman & Sen (2009):
    LOD = (n/2) * log10(1 + F * df_reg / df_res)
    
    LOD is always >= 0 since F >= 0.
    """
    if df_res <= 0 or np.isnan(f_stat) or f_stat < 0:
        return 0.0  # Return 0 for invalid cases
    return (n / 2) * np.log10(1 + f_stat * df_reg / df_res)


def f_to_lod_vectorized(f_stats, n, df_reg, df_res):
    """Vectorized F to LOD conversion. LOD is always >= 0."""
    if df_res <= 0:
        return np.zeros_like(f_stats)
    
    # Ensure F >= 0
    f_stats = np.maximum(f_stats, 0)
    
    result = (n / 2) * np.log10(1 + f_stats * df_reg / df_res)
    
    # Handle any remaining NaN/inf
    result = np.nan_to_num(result, nan=0.0, posinf=100.0, neginf=0.0)
    
    return result


def lod_to_pve(lod, n):
    """
    Estimate PVE from LOD score.
    Formula: PVE = 1 - 10^(-2*LOD/n)
    """
    if np.isnan(lod) or n <= 0 or lod <= 0:
        return 0.0
    pve = 1 - 10**(-2 * lod / n)
    return max(0, min(1, pve))


def compute_block_lod_permutations(X, y, n_perms, rng):
    """
    Fast vectorized permutation test for a single block using hat matrix.
    
    Returns observed LOD and LOD for all permutations.
    This is the key to speed - vectorize all permutations at once per block.
    
    Parameters:
    -----------
    X : array (n_samples, n_founders)
        Design matrix for this block
    y : array (n_samples,)
        Phenotype values
    n_perms : int
        Number of permutations
    rng : numpy.random.Generator
        Random number generator
        
    Returns:
    --------
    obs_f : float
        Observed F-statistic
    obs_lod : float
        Observed LOD score
    perm_lods : array (n_perms,)
        LOD scores from all permutations
    """
    n, p = X.shape

    # Add intercept
    X1 = np.column_stack([np.ones(n), X])

    # Degrees of freedom
    df_reg = p
    df_res = n - p - 1

    if df_res <= 0:
        return 0.0, 0.0, np.zeros(n_perms)

    # Precompute projection (hat) matrix: H = X @ (X'X)^-1 @ X'
    try:
        XtX = X1.T @ X1
        XtX_inv = np.linalg.inv(XtX)
        H = X1 @ XtX_inv @ X1.T
    except np.linalg.LinAlgError:
        return 0.0, 0.0, np.zeros(n_perms)

    # Total sum of squares (constant across permutations)
    y_mean = y.mean()
    ss_total = np.sum((y - y_mean) ** 2)

    if ss_total <= 0:
        return 0.0, 0.0, np.zeros(n_perms)

    # Observed model
    y_hat = H @ y
    ss_res = np.sum((y - y_hat) ** 2)
    ss_reg = ss_total - ss_res
    
    # Clamp ss_reg to be non-negative (numerical precision fix)
    ss_reg = max(ss_reg, 0)

    if ss_res <= 0:
        obs_f = np.inf
    else:
        obs_f = (ss_reg / df_reg) / (ss_res / df_res)

    obs_lod = f_to_lod(obs_f, n, df_reg, df_res)

    # ---- Vectorized permutations ----
    # Generate all permutation indices at once
    perm_indices = np.array([rng.permutation(n) for _ in range(n_perms)])
    Yp = y[perm_indices]  # (n_perms, n)

    # Compute predicted values for all permutations at once
    Yp_hat = Yp @ H.T  # (n_perms, n)

    # Residual sum of squares for all permutations
    ss_res_perm = np.sum((Yp - Yp_hat) ** 2, axis=1)
    
    # Clamp ss_res to be at most ss_total (numerical precision)
    ss_res_perm = np.minimum(ss_res_perm, ss_total)
    
    ss_reg_perm = ss_total - ss_res_perm

    # F-statistics for all permutations
    # Clamp ss_reg to be non-negative (numerical precision fix)
    ss_reg_perm = np.maximum(ss_reg_perm, 0)
    
    with np.errstate(divide='ignore', invalid='ignore'):
        f_perm = (ss_reg_perm / df_reg) / (ss_res_perm / df_res)
    
    # Clamp F to be non-negative (safety)
    f_perm = np.maximum(f_perm, 0)

    # Convert to LOD
    perm_lods = f_to_lod_vectorized(f_perm, n, df_reg, df_res)

    return obs_f, obs_lod, perm_lods


def test_phenotype(pheno, merged, blocks_dict, block_cols, n_perms=1000):
    """
    Test all blocks for a single phenotype.
    
    Strategy for speed + correctness:
    1. For each block, compute observed LOD and all n_perms permutation LODs (vectorized)
    2. For each permutation i, take max LOD across all blocks -> genome-wide null distribution
    3. Threshold = 95th percentile of max LOD distribution
    """
    print(f"\nTesting phenotype: {pheno}")

    # Get phenotype values (remove NaN)
    pheno_data = merged[[pheno] + block_cols].dropna(subset=[pheno])

    if len(pheno_data) < 10:
        print(f"  Skipping {pheno}: too few samples ({len(pheno_data)})")
        return [], {}, np.array([])

    y = pheno_data[pheno].values
    n = len(y)

    print(f"  N samples: {n}")
    print(f"  Testing {len(blocks_dict)} blocks with {n_perms} permutations...")

    # Build X matrices for all valid blocks
    X_matrices = []
    block_ids = []
    block_infos = []

    for block_id, block_info in blocks_dict.items():
        founder_cols = block_info['columns']
        X = pheno_data[founder_cols].values

        # Skip blocks with no variation
        if X.sum() == 0 or np.all(X == X[0]):
            continue

        X_matrices.append(X)
        block_ids.append(block_id)
        block_infos.append(block_info)

    n_blocks = len(X_matrices)
    if n_blocks == 0:
        print(f"  No valid blocks for {pheno}")
        return [], {}, np.array([])

    print(f"  Valid blocks: {n_blocks}")

    # Use same RNG for all blocks so permutation i uses same y permutation across blocks
    rng = np.random.default_rng(12345)

    # Store results
    obs_f_all = np.zeros(n_blocks)
    obs_lod_all = np.zeros(n_blocks)
    perm_lods_all = np.zeros((n_blocks, n_perms))  # (n_blocks, n_perms)

    print(f"  Running vectorized permutation tests...")

    # For each block, compute observed and permutation LODs
    for i, X in enumerate(X_matrices):
        # IMPORTANT: Reset RNG state for each block so permutation j 
        # corresponds to the same phenotype permutation across all blocks
        block_rng = np.random.default_rng(12345)
        
        obs_f, obs_lod, perm_lods = compute_block_lod_permutations(
            X, y, n_perms, block_rng
        )
        obs_f_all[i] = obs_f
        obs_lod_all[i] = obs_lod
        perm_lods_all[i, :] = perm_lods

    # For each permutation, get the MAX LOD across all blocks
    # This gives the genome-wide null distribution (Churchill & Doerge 1994)
    max_lod_per_perm = np.max(perm_lods_all, axis=0)  # (n_perms,)

    # Calculate genome-wide thresholds
    thresholds = {
        'LOD_thresh_0.05': np.percentile(max_lod_per_perm, 95),
        'LOD_thresh_0.01': np.percentile(max_lod_per_perm, 99),
        'LOD_thresh_0.10': np.percentile(max_lod_per_perm, 90),
        'LOD_thresh_0.001': np.percentile(max_lod_per_perm, 99.9)
    }

    print(f"  Thresholds: LOD > {thresholds['LOD_thresh_0.05']:.2f} (5% FWER), "
          f"LOD > {thresholds['LOD_thresh_0.01']:.2f} (1% FWER)")

    # Compile results
    results = []
    for i, (block_id, block_info) in enumerate(zip(block_ids, block_infos)):
        obs_f = obs_f_all[i]
        obs_lod = obs_lod_all[i]
        pve = lod_to_pve(obs_lod, n)

        # Significance based on genome-wide permutation thresholds
        sig_0_05 = obs_lod > thresholds['LOD_thresh_0.05'] if not np.isnan(obs_lod) else False
        sig_0_01 = obs_lod > thresholds['LOD_thresh_0.01'] if not np.isnan(obs_lod) else False

        results.append({
            'chrom': block_info['chrom'],
            'start_pos': block_info['start_pos'],
            'end_pos': block_info['end_pos'],
            'block_id': block_id,
            'phenotype': pheno,
            'n_samples': n,
            'n_founders': len(block_info['columns']),
            'f_stat': obs_f,
            'lod': obs_lod,
            'pve': pve,
            'sig_0.05': sig_0_05,
            'sig_0.01': sig_0_01
        })

    # Summary
    n_sig_05 = sum(1 for r in results if r['sig_0.05'])
    n_sig_01 = sum(1 for r in results if r['sig_0.01'])
    print(f"  Completed: {len(results)} blocks tested")
    print(f"  Significant at 5% FWER: {n_sig_05}")
    print(f"  Significant at 1% FWER: {n_sig_01}")

    return results, thresholds, max_lod_per_perm


def main():
    parser = argparse.ArgumentParser(
        description='QTL mapping with permutation testing - Output: F-stat, LOD, PVE'
    )
    parser.add_argument('--pheno-start', type=int, required=True,
                        help='Starting phenotype index')
    parser.add_argument('--pheno-end', type=int, required=True,
                        help='Ending phenotype index')
    parser.add_argument('--n-perms', type=int, default=1000,
                        help='Number of permutations for threshold calculation')
    parser.add_argument('--output-dir', type=str, default='qtl_results_parts',
                        help='Output directory')

    args = parser.parse_args()

    print("=" * 80)
    print("QTL Mapping with Permutation Testing (FAST + CORRECT)")
    print("Following Churchill & Doerge (1994) for genome-wide threshold")
    print("Output: F-stat, LOD, PVE (effect size)")
    print(f"Phenotypes: {args.pheno_start} to {args.pheno_end}")
    print(f"Permutations: {args.n_perms}")
    print("=" * 80)

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    print("\nReading phenotype data...")
    pheno_df = pd.read_csv('FemaleHeadExpression.txt', sep='\t')
    pheno_df['id'] = pheno_df['patRIL']

    # Read filtered phenotype list
    print("Reading filtered phenotype list...")
    if not os.path.exists('phenotypes_to_test.csv'):
        print("ERROR: phenotypes_to_test.csv not found!")
        print("Run generate_slurm_jobs.py first to create filtered phenotype list.")
        sys.exit(1)

    pheno_list_df = pd.read_csv('phenotypes_to_test.csv')
    pheno_cols = pheno_list_df['phenotype'].tolist()
    print(f"Loaded {len(pheno_cols)} phenotypes to test (filtered by h2)")

    # Verify these phenotypes exist in the expression data
    pheno_cols = [p for p in pheno_cols if p in pheno_df.columns]
    print(f"Phenotypes available in expression data: {len(pheno_cols)}")

    print("\nReading block matrix...")
    block_df = pd.read_csv('Cross_RILS_Hap_Blocks.csv', low_memory=False)

    print("Merging datasets...")
    merged = block_df.merge(pheno_df, on='id', how='inner')
    print(f"Merged data shape: {merged.shape}")

    # Select subset of phenotypes for this job
    pheno_subset = pheno_cols[args.pheno_start:args.pheno_end]
    print(f"Testing phenotypes: {pheno_subset[:5]}... ({len(pheno_subset)} total)")

    # Group columns by block
    print("\nGrouping founder columns by block...")
    block_cols = [col for col in block_df.columns if col != 'id']
    blocks_dict = {}

    for col in block_cols:
        match = re.match(r'(.+?)_(\d+)(?:-(\d+))?_(AA|BB)(\d+)', col)
        if match:
            chrom = match.group(1)
            start_pos = int(match.group(2))
            end_pos = int(match.group(3)) if match.group(3) else start_pos

            block_id = f"{chrom}_{start_pos}-{end_pos}" if start_pos != end_pos else f"{chrom}_{start_pos}"

            if block_id not in blocks_dict:
                blocks_dict[block_id] = {
                    'chrom': chrom,
                    'start_pos': start_pos,
                    'end_pos': end_pos,
                    'columns': []
                }
            blocks_dict[block_id]['columns'].append(col)

    print(f"Found {len(blocks_dict)} unique blocks")

    # Test each phenotype
    all_results = []
    all_thresholds = {}

    for pheno in pheno_subset:
        pheno_results, pheno_thresholds, max_lod_dist = test_phenotype(
            pheno, merged, blocks_dict, block_cols, n_perms=args.n_perms
        )
        all_results.extend(pheno_results)

        if pheno_thresholds:
            all_thresholds[pheno] = pheno_thresholds

            # Save permutation distribution
            if len(max_lod_dist) > 0:
                perm_file = os.path.join(args.output_dir, f'perm_max_lod_{pheno}.npy')
                np.save(perm_file, max_lod_dist)

    # Save results
    if len(all_results) > 0:
        results_df = pd.DataFrame(all_results)
        results_df = results_df.sort_values('lod', ascending=False)

        output_file = os.path.join(
            args.output_dir,
            f'qtl_results_{args.pheno_start}_{args.pheno_end}.csv'
        )
        results_df.to_csv(output_file, index=False)

        print(f"\n{'=' * 80}")
        print(f"Saved results to {output_file}")
        print(f"Total tests: {len(results_df)}")
        print(f"Significant at 5% FWER: {results_df['sig_0.05'].sum()}")
        print(f"Significant at 1% FWER: {results_df['sig_0.01'].sum()}")

        print(f"\nLOD score summary:")
        print(f"  Max LOD: {results_df['lod'].max():.2f}")
        print(f"  Mean LOD: {results_df['lod'].mean():.2f}")
        print(f"  Median LOD: {results_df['lod'].median():.2f}")

        print(f"\nTop hits (significant at 5% FWER):")
        sig_hits = results_df[results_df['sig_0.05']].head(10)
        for _, row in sig_hits.iterrows():
            print(f"  {row['block_id']} ({row['phenotype']}): "
                  f"LOD={row['lod']:.2f}, PVE={row['pve'] * 100:.1f}%")

        print(f"{'=' * 80}")

        # Save thresholds
        if all_thresholds:
            thresh_file = os.path.join(
                args.output_dir,
                f'thresholds_{args.pheno_start}_{args.pheno_end}.csv'
            )
            thresh_df = pd.DataFrame(all_thresholds).T
            thresh_df.index.name = 'phenotype'
            thresh_df.to_csv(thresh_file)
            print(f"Saved thresholds to {thresh_file}")
    else:
        print("No results generated!")


if __name__ == '__main__':
    main()
