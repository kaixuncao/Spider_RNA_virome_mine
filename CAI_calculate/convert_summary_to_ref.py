import pandas as pd
import os

def convert_summary_to_cai_refs(input_csv, output_dir="reference_tables"):
    """
    Convert the summary table of codon counts into a single reference table (.txt),
    For use with calc_cai_coords2.py.
    """
    
    # 1. Read CSV
    print(f"Reading {input_csv}...")
    try:
        df = pd.read_csv(input_csv)
    except FileNotFoundError:
        print("Error: The input file cannot be found. Please ensure that the CSV file is in the current directory.")
        return

    # 2. Create an output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    # 3. Define 64 codons (ensuring that the column names match)
    bases = ['T', 'C', 'A', 'G']
    codons = [a+b+c for a in bases for b in bases for c in bases]
    
    # Check whether the CSV file contains all the codon columns
    missing_cols = [c for c in codons if c not in df.columns]
    if missing_cols:
        print(f"Error: The CSV file is missing the following codon columns: {missing_cols}")
        return

    # 4. Iterate through each row (each species/sample)
    count = 0
    for index, row in df.iterrows():
        # Retrieve the sample ID (to be used as the file name)
        # Assume the ID column is named 'ID'; if not, attempt to identify it automatically
        sample_id = str(row.get('Filename', f"Sample_{index}"))
        
        # Remove invalid characters (such as colons) from file names
        safe_filename = sample_id.replace(":", "_").replace("/", "_").replace("\\", "_")
        output_file = os.path.join(output_dir, f"{safe_filename}.txt")
        
        # Write to a file
        with open(output_file, 'w') as f:
            # Include a header (optional, but recommended for clarity)
            f.write("Codon\tCount\n")
            
            for codon in codons:
                # Obtain the count for this codon
                codon_count = row[codon]
                # Write: TTT <tab> 123
                f.write(f"{codon}\t{codon_count}\n")
        
        count += 1

    print(f"Conversion complete!")
    print(f"Successfully generated {count} reference files in '{output_dir}/'.")

if __name__ == "__main__":
    # Replace the filename here with the actual name of the CSV file you have generated
    input_csv_name = "all_Final_Codon_Usage_Summary.csv" 
    
    convert_summary_to_cai_refs(input_csv_name)
