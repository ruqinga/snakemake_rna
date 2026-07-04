import pandas as pd
from pathlib import Path

# config: when called by Snakemake
if 'snakemake' in globals():
    input_files = snakemake.input  # Input files from Snakemake
    output_path = Path(snakemake.output.merged_fpkm)  # Output path from Snakemake
else:
    # Direct Python execution: setting the paths manually
    input_files = list(Path('./Results/06_norm').rglob('genes.fpkm_tracking'))
    output_path = Path('./Results/06_norm/merged_fpkm.txt')

def merge_fpkm_files():
    merged_df = None

    for fp in sorted(input_files):
        # Safely read the file, skipping comment lines
        df = pd.read_csv(fp, sep='\t', comment='#')

        # Extract the sample ID from the folder name
        sample_id = Path(fp).parent.name  # Get the parent folder name
        sample_id = sample_id.replace('_pe', '').replace('_se', '')  # Remove _pe/_se suffixes

        # Ensure the necessary columns exist
        required_columns = ['gene_short_name', 'FPKM']
        if not all(col in df.columns for col in required_columns):
            raise ValueError(f"File {fp} is missing required columns: {required_columns}")

        # Extract necessary columns and rename the expression column
        df = df[required_columns].copy()
        df.columns = ['gene', sample_id]

        # Average gene entries if duplicates exist
        df = df.groupby('gene').mean().reset_index()

        # Merge data
        if merged_df is None:
            merged_df = df
        else:
            merged_df = pd.merge(
                merged_df,
                df,
                on='gene',
                how='outer',
                validate='one_to_one'
            )

    # Sort and save the final result
    merged_df.sort_values('gene').to_csv(
        output_path,
        sep='\t',
        index=False,
        encoding='utf-8'
    )

if __name__ == '__main__':
    # Create the output directory if it doesn't exist
    output_path.parent.mkdir(parents=True, exist_ok=True)
    merge_fpkm_files()
