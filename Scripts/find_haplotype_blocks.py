import pandas as pd
import numpy as np
import re
from collections import defaultdict
import sys

print("Reading Cross_RILS_Hap_Probs.csv with pandas...")
# Read with pandas - more forgiving of parsing issues
df = pd.read_csv('Cross_RILS_Hap_Probs.csv', low_memory=False)

print(f"Data shape: {df.shape}")
print(f"Number of individuals: {len(df)}")
print(f"Number of columns: {len(df.columns)}")

# Check for any NaN in id column
if df['id'].isna().any():
    print("WARNING: Found NaN values in id column")
    df = df.dropna(subset=['id'])
    print(f"After dropping NaN ids: {len(df)} rows")

# Parse column names to extract positions
print("\nParsing positions from column names...")
position_info = {}  # {position_key: {'chrom': ..., 'pos': ..., 'aa_cols': [...], 'bb_cols': [...]}}

for col in df.columns:
    if col == 'id':
        continue
    
    # Match pattern: chromosome_position_founder (e.g., 2L_100000_AA1)
    match = re.match(r'(.+?)_(\d+)_(AA|BB)(\d+)', col)
    if match:
        chrom = match.group(1)
        pos = int(match.group(2))
        founder_type = match.group(3)
        
        pos_key = f"{chrom}_{pos}"
        
        if pos_key not in position_info:
            position_info[pos_key] = {
                'chrom': chrom,
                'pos': pos,
                'aa_cols': [],
                'bb_cols': []
            }
        
        if founder_type == 'AA':
            position_info[pos_key]['aa_cols'].append(col)
        else:
            position_info[pos_key]['bb_cols'].append(col)

print(f"Found {len(position_info)} positions")

# Sort positions by chromosome and position
sorted_positions = sorted(position_info.items(), 
                         key=lambda x: (x[1]['chrom'], x[1]['pos']))

print("\nFinding likely founders at each position for each individual...")

position_signatures = []

for idx, (pos_key, info) in enumerate(sorted_positions):
    if idx % 1000 == 0:
        print(f"Processing position {idx}/{len(sorted_positions)}...")
    
    # Get AA and BB columns for this position
    aa_cols = info['aa_cols']
    bb_cols = info['bb_cols']
    
    if not aa_cols or not bb_cols:
        continue
    
    # For each individual, find which AA and BB founder(s) have max probability
    aa_data = df[aa_cols].values
    bb_data = df[bb_cols].values
    
    # Store the signature for each individual at this position
    individual_signatures = []
    
    for i in range(len(df)):
        aa_row = aa_data[i]
        bb_row = bb_data[i]
        
        aa_max_val = np.max(aa_row)
        bb_max_val = np.max(bb_row)
        
        # Find which AA founders have max value (handles ties)
        likely_aa = tuple(sorted([j+1 for j, val in enumerate(aa_row) if val == aa_max_val]))
        likely_bb = tuple(sorted([j+1 for j, val in enumerate(bb_row) if val == bb_max_val]))
        
        individual_signatures.append((likely_aa, likely_bb))
    
    position_signatures.append({
        'chrom': info['chrom'],
        'pos': info['pos'],
        'pos_key': pos_key,
        'signatures': tuple(individual_signatures)  # One signature per individual
    })

print(f"Analyzed {len(position_signatures)} positions")

# Find consecutive regions where ALL individuals maintain their own consistent founders
print("\nIdentifying haplotype blocks...")

regions = []
current_region = None

for i, pos_info in enumerate(position_signatures):
    if current_region is None:
        # Start new region
        current_region = {
            'chrom': pos_info['chrom'],
            'start_pos': pos_info['pos'],
            'end_pos': pos_info['pos'],
            'num_pos': 1,
            'signatures': pos_info['signatures']
        }
    elif (pos_info['chrom'] == current_region['chrom'] and 
          pos_info['signatures'] == current_region['signatures']):
        # Same chromosome AND all individuals have same founders as previous position
        # Extend current region
        current_region['end_pos'] = pos_info['pos']
        current_region['num_pos'] += 1
    else:
        # Different chromosome or at least one individual changed founders
        # Save current region and start new one
        regions.append(current_region)
        current_region = {
            'chrom': pos_info['chrom'],
            'start_pos': pos_info['pos'],
            'end_pos': pos_info['pos'],
            'num_pos': 1,
            'signatures': pos_info['signatures']
        }

# Don't forget the last region
if current_region is not None:
    regions.append(current_region)

print(f"Found {len(regions)} haplotype blocks")

# Create output dataframe
output_df = pd.DataFrame({
    'chrom': [r['chrom'] for r in regions],
    'start_pos': [r['start_pos'] for r in regions],
    'end_pos': [r['end_pos'] for r in regions],
    'num_pos': [r['num_pos'] for r in regions]
})

# Save to CSV
output_file = 'haplotype_blocks.csv'
output_df.to_csv(output_file, index=False)

print(f"\nSaved results to {output_file}")
print(f"\nFirst 10 regions:")
print(output_df.head(10))

print(f"\nSummary by chromosome:")
summary = output_df.groupby('chrom').agg({
    'num_pos': ['count', 'sum', 'mean', 'max']
})
summary.columns = ['num_blocks', 'total_positions', 'avg_block_size', 'max_block_size']
print(summary)

print(f"\nOverall statistics:")
print(f"Total blocks: {len(regions)}")
print(f"Total positions analyzed: {len(position_signatures)}")
print(f"Average block size: {output_df['num_pos'].mean():.2f} positions")
print(f"Median block size: {output_df['num_pos'].median():.0f} positions")
print(f"Max block size: {output_df['num_pos'].max()} positions")
