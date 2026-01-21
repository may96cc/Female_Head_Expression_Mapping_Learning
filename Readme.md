# Female Head Expression

QTL mapping of gene expression traits in female *Drosophila melanogaster* heads using the Drosophila Synthetic Population Resource (DSPR). This project identifies genetic loci associated with variation in gene expression using recombinant inbred line (RIL) crosses. Modeling via deep learning (MLPs) to be added

Data from original publication : King, E. G., Sanderson, B. J., McNeil, C. L., Long, A. D., & Macdonald, S. J. (2014). Genetic dissection of the Drosophila melanogaster female head transcriptome reveals widespread allelic heterogeneity. PLoS Genetics, 10(5), e1004322.

DSPR website : https://wfitch.bio.uci.edu/~dspr/

## Project Structure

```
Female_Head_Expression/
├── RawData/
│   ├── FemaleHeadExpression.txt        # Gene expression measurements
│   ├── h2.csv                          # Heritability estimates
│   └── cross_hap_phenotype.csv         # Cross haplotype phenotype data
├── ProcessedData/
│   ├── A_RILS_Hap_Probs.csv            # A-panel RIL haplotype probabilities
│   ├── B_RILS_Hap_Probs.csv            # B-panel RIL haplotype probabilities
│   ├── Cross_RILS_Hap_Blocks.csv       # Cross RIL haplotype blocks
│   └── Cross_RILS_Hap_Probs.csv        # Cross RIL haplotype probabilities
├── Scripts/
│   ├── make_hap_prob_mat_A.R           # Generate A-panel haplotype matrix
│   ├── make_hap_prob_mat_A.sh          # SLURM submission script
│   ├── make_hap_prob_mat_B.R           # Generate B-panel haplotype matrix
│   ├── make_hap_prob_mat_B.sh          # SLURM submission script
│   ├── find_haplotype_blocks.py        # Identify haplotype blocks
│   ├── find_haplotype_blocks.sh        # SLURM submission script
│   └── qtl_mapping_fast.py             # QTL mapping with permutation testing
└── Readme.md
```

## Data Description

### RawData/

#### FemaleHeadExpression.txt
Tab-delimited file containing gene expression values for RIL crosses.

| Column | Description |
|--------|-------------|
| patRIL | Paternal RIL identifier |
| matRIL | Maternal RIL identifier |
| CG##### | Expression values for each gene (one column per gene) |

#### h2.csv
Comma-separated file containing heritability (h²) estimates for each gene.

| Column | Description |
|--------|-------------|
| (row names) | Gene identifier (e.g., CG17210) |
| x | Heritability estimate |

#### cross_hap_phenotype.csv
Comma-separated file containing cross and haplotype phenotype information.

| Column | Description |
|--------|-------------|
| patRIL | Paternal RIL identifier |
| matRIL | Maternal RIL identifier |
| sex | Sex of the sample (F = Female) |
| id | Cross identifier |

### ProcessedData/

#### A_RILS_Hap_Probs.csv
Haplotype probabilities for A-panel RILs.

| Column | Description |
|--------|-------------|
| ril | RIL identifier |
| {chr}_{pos}_{hap} | Haplotype probability (e.g., `2L_100000_AA1` = chromosome 2L, position 100000, founder AA1) |

#### B_RILS_Hap_Probs.csv
Haplotype probabilities for B-panel RILs.

| Column | Description |
|--------|-------------|
| ril | RIL identifier |
| {chr}_{pos}_{hap} | Haplotype probability (e.g., `2L_100000_BB1` = chromosome 2L, position 100000, founder BB1) |

#### Cross_RILS_Hap_Blocks.csv
Haplotype block assignments for cross RILs. Values represent discrete block calls.

| Column | Description |
|--------|-------------|
| id | Cross RIL identifier |
| {chr}_{pos}_{hap} | Haplotype block assignment (binary/discrete values) |

#### Cross_RILS_Hap_Probs.csv
Haplotype probabilities for cross RILs.

| Column | Description |
|--------|-------------|
| id | Cross RIL identifier |
| {chr}_{pos}_{hap} | Haplotype probability at each genomic position |

### Scripts/

#### make_hap_prob_mat_A.R
Extracts haplotype probabilities for A-panel RILs from the `DSPRqtlDataA` R package. Merges all genomic positions into a single matrix with columns named `{chrom}_{position}_AA{1-8}`.

**Output:** `A_RILS_Hap_Probs.csv`

#### make_hap_prob_mat_B.R
Extracts haplotype probabilities for B-panel RILs from the `DSPRqtlDataB` R package. Merges all genomic positions into a single matrix with columns named `{chrom}_{position}_BB{1-8}`.

**Output:** `B_RILS_Hap_Probs.csv`

#### find_haplotype_blocks.py
Identifies haplotype blocks from cross RIL haplotype probabilities. For each genomic position, determines the most likely founder haplotype for each individual, then groups consecutive positions with identical founder signatures into blocks.

**Input:** `Cross_RILS_Hap_Probs.csv`  
**Output:** `haplotype_blocks.csv`

#### qtl_mapping_fast.py
Performs QTL mapping with permutation testing following Churchill & Doerge (1994). Uses vectorized permutation tests with the hat matrix trick for computational efficiency.

**Features:**
- Computes F-statistics, LOD scores, and PVE (percent variance explained)
- Genome-wide significance thresholds via permutation (5%, 1%, 0.1% FWER)
- Parallelizable across phenotypes via command-line arguments

**Usage:**
```bash
python qtl_mapping_fast.py --pheno-start 0 --pheno-end 100 --n-perms 1000 --output-dir qtl_results_parts
```

**Inputs:** `FemaleHeadExpression.txt`, `Cross_RILS_Hap_Blocks.csv`, `phenotypes_to_test.csv`  
**Outputs:** `qtl_results_{start}_{end}.csv`, `thresholds_{start}_{end}.csv`

#### Shell scripts (*.sh)
SLURM batch submission scripts for running the corresponding R/Python scripts on the cluster.

## Methods

### DSPR Overview
The Drosophila Synthetic Population Resource (DSPR) consists of two panels of RILs (A and B), each derived from 8 founder lines. This project uses crosses between the A and B panels to map expression QTLs (eQTLs).

### Analysis Pipeline
1. **Haplotype probability extraction** - Extract founder haplotype probabilities from DSPR R packages (`DSPRqtlDataA`, `DSPRqtlDataB`)
2. **Haplotype block identification** - Group consecutive genomic positions with consistent founder assignments
3. **QTL mapping** - Test association between haplotype blocks and gene expression using F-tests
4. **Significance thresholds** - Genome-wide thresholds determined by permutation testing (Churchill & Doerge 1994)

### Statistical Framework
- **Test statistic:** F-statistic from regression of expression on founder haplotype probabilities
- **Effect size:** LOD score and PVE (percent variance explained)
- **Multiple testing correction:** Permutation-based family-wise error rate (FWER) control

## Usage

### Dependencies
- **R packages:** `DSPRqtlDataA`, `DSPRqtlDataB`, `data.table`
- **Python packages:** `pandas`, `numpy`, `scipy`

### Running the Pipeline

1. Generate haplotype probability matrices:
```bash
sbatch Scripts/make_hap_prob_mat_A.sh
sbatch Scripts/make_hap_prob_mat_B.sh
```

2. Identify haplotype blocks:
```bash
sbatch Scripts/find_haplotype_blocks.sh
```

3. Run QTL mapping (parallelized by phenotype batches):
```bash
python Scripts/qtl_mapping_fast.py --pheno-start 0 --pheno-end 100 --n-perms 1000
```

## Contact

[Your name/email]

## References

- King EG, Macdonald SJ, Long AD (2012). Properties and power of the Drosophila Synthetic Population Resource for the routine dissection of complex traits. *Genetics* 191(3):935-949.
- Churchill GA, Doerge RW (1994). Empirical threshold values for quantitative trait mapping. *Genetics* 138(3):963-971.
- Broman KW, Sen S (2009). *A Guide to QTL Mapping with R/qtl*. Springer.
