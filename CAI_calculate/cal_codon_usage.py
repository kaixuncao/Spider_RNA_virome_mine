import os
import pandas as pd
from Bio import SeqIO, Entrez
from collections import Counter

# Configuration section: please use your NCBI email address and API

Entrez.email = "NCBI email address" 
Entrez.api_key = "API KEY" 

# ORF prediction threshold: retain only those longer than 300
MIN_ORF_LENGTH = 300 

# Global cache
TAXON_CACHE = {}

# Functional Modules: APIs and Computing

def fetch_genetic_code(tax_id):
    """Retrieve the genetic code table ID for a species via the NCBI API"""
    if not tax_id: return 1
    if tax_id in TAXON_CACHE: return TAXON_CACHE[tax_id]
    
    try:
        print(f"  [API] Querying TaxID {tax_id}...")
        handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        if records:
            # Obtain the nuclear genome codon table (GCId)
            gc_id = int(records[0]['GeneticCode']['GCId'])
            TAXON_CACHE[tax_id] = gc_id
            return gc_id
    except Exception:
        pass
    return 1 # On failure, the standard error code is returned by default

def calculate_gc3s(seq):
    """calculate GC3s"""
    clean_seq = seq[:len(seq)-(len(seq)%3)]
    if not clean_seq: return 0
    third_pos = clean_seq[2::3]
    gc = third_pos.count('G') + third_pos.count('C') + third_pos.count('g') + third_pos.count('c')
    return (gc / len(third_pos) * 100) if len(third_pos) > 0 else 0

def get_codon_counts(seq):
    """Retrieve the number of codons"""
    cnt = Counter()
    clean_seq = str(seq[:len(seq)-(len(seq)%3)]).upper()
    codons = [clean_seq[i:i+3] for i in range(0, len(clean_seq), 3)]
    cnt.update(codons)
    return cnt

# Core module: ORF prediction (for unannotated files)

def find_orfs_in_sequence(dna_seq, table_id, min_len):
    """
    Search for potential ORFs (open reading frames) in DNA sequences.
    Logic: Translate the three forward reading frames, looking for an M...* structure.
    """
    orfs = []
    seq_len = len(dna_seq)
    
    # Iterate through the three forward decoding boxes (0, 1, 2)
    # Reverse-complementary sequences are not currently taken into account; WGS contigs are typically submitted in the forward orientation, and in order to minimise noise
    for frame in range(3):
        # Extract the relevant sections for translation
        coding_dna = dna_seq[frame:]
        coding_dna = coding_dna[:len(coding_dna) - (len(coding_dna) % 3)]
        
        if len(coding_dna) < min_len: continue
        
        # Translation
        try:
            protein = coding_dna.translate(table=table_id)
        except Exception:
            protein = coding_dna.translate(table=1) # Return to the standard table

        # Searching for M (Start) and * (Stop) in protein sequences
        prot_str = str(protein)
        aa_len = len(prot_str)
        
        start_idx = -1
        
        for i, aa in enumerate(prot_str):
            if aa == 'M' and start_idx == -1:
                start_idx = i
            elif aa == '*':
                if start_idx != -1:
                    # Find a complete ORF
                    orf_aa_len = i - start_idx
                    orf_dna_len = orf_aa_len * 3
                    
                    if orf_dna_len >= min_len:
                        # Convert back to DNA coordinates and extract the sequence
                        # DNA start = frame + (start_idx * 3)
                        # DNA end = DNA start + orf_dna_len
                        dna_start = frame + (start_idx * 3)
                        dna_end = dna_start + orf_dna_len + 3 # +3 Contains a stop codon
                        
                        extracted_seq = dna_seq[dna_start:dna_end]
                        orfs.append(extracted_seq)
                    
                    # Reset `start_idx` and search for the next ORF
                    start_idx = -1 
                # If * is encountered and there is no start_idx, ignore it
    
    return orfs


# Main Process


def process_files_force_mode(root_dir):
    results = []
    
    for dirpath, dirnames, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.endswith(".gbff") or filename.endswith(".gbk"):
                full_path = os.path.join(dirpath, filename)
                parent_folder = os.path.basename(dirpath)
                file_stem = os.path.splitext(filename)[0]
                sample_id = f"{parent_folder}:{file_stem}"
                
                print(f"Processing: {filename} ...")
                
                final_cds_list = []
                extraction_method = "Unknown"
                genetic_code = 1
                tax_id = None
                
                try:
                    # Retrieve all records
                    records = list(SeqIO.parse(full_path, "genbank"))
                    if not records: continue
                    
                    # 1. Attempt to retrieve the TaxID and genetic code
                    for feat in records[0].features:
                        if feat.type == "source":
                            if 'db_xref' in feat.qualifiers:
                                for xref in feat.qualifiers['db_xref']:
                                    if xref.startswith("taxon:"):
                                        tax_id = xref.split(":")[1]
                    
                    # Verifying the genetic code online
                    genetic_code = fetch_genetic_code(tax_id) if tax_id else 1
                    
                    # 2. Strategy A: Attempt a standard CDS extraction
                    for record in records:
                        for feat in record.features:
                            if feat.type == "CDS":
                                try:
                                    seq = feat.location.extract(record.seq)
                                    if len(seq) % 3 == 0 and len(seq) > 0:
                                        final_cds_list.append(seq)
                                except: pass
                    
                    if final_cds_list:
                        extraction_method = "Annotation_CDS"
                        print(f"  -> Found {len(final_cds_list)} CDS via Annotation.")
                    
                    # 3. Strategy B: If no CDS is available, initiate ORF prediction (aggressive mode)
                    else:
                        print(f"  -> No Annotation found. Scanning for ORFs (Min {MIN_ORF_LENGTH}bp)...")
                        extraction_method = "Predicted_ORF"
                        
                        for record in records:
                            # Retrieve the sequence (convert to a Seq object)
                            dna_seq = record.seq
                            # Predict
                            orfs = find_orfs_in_sequence(dna_seq, genetic_code, MIN_ORF_LENGTH)
                            final_cds_list.extend(orfs)
                        
                        print(f"  -> Predicted {len(final_cds_list)} putative ORFs.")

                    # 4. If you can’t even find a predict, skip this
                    if not final_cds_list:
                        print("  -> Failed to extract valid sequences.")
                        continue
                        
                    # 5. Save FASTA and calculate metrics
                    out_fasta = os.path.join(dirpath, f"{file_stem}_CDS.fasta")
                    all_counts = Counter()
                    gc3s_sum = 0
                    
                    with open(out_fasta, "w") as f:
                        for i, seq in enumerate(final_cds_list):
                            # ID format: Sample name_Method_Serial number
                            seq_id = f"{sample_id}_{extraction_method}_{i+1}"
                            f.write(f">{seq_id}\n{str(seq)}\n")
                            
                            all_counts.update(get_codon_counts(seq))
                            gc3s_sum += calculate_gc3s(seq)
                    
                    avg_gc3s = gc3s_sum / len(final_cds_list)
                    
                    # 6. Write summary data
                    row = {
                        "ID": sample_id,
                        "Folder": parent_folder,
                        "Filename": filename,
                        "Method": extraction_method, # Determine whether it is a real CDS or a forecast
                        "Genetic_Code": genetic_code,
                        "Total_Seqs": len(final_cds_list),
                        "Avg_GC3s": round(avg_gc3s, 2)
                    }
                    # Add 64 codons usage
                    bases = ['T', 'C', 'A', 'G']
                    codons = [a+b+c for a in bases for b in bases for c in bases]
                    for codon in codons:
                        row[codon] = all_counts.get(codon, 0)
                        
                    results.append(row)

                except Exception as e:
                    print(f"  Error: {e}")
                    
    return results

if __name__ == "__main__":
    print("Starting Force-Extraction Mode...")
    data = process_files_force_mode(".")
    
    if data:
        df = pd.DataFrame(data)
        # Sort the columns
        first_cols = ['ID', 'Folder', 'Filename', 'Method', 'Genetic_Code', 'Total_Seqs', 'Avg_GC3s']
        cols = first_cols + [c for c in df.columns if c not in first_cols]
        df = df[cols]
        
        df.to_csv("Final_Codon_Usage_Summary.csv", index=False)
        print("\nAll Done! Check 'Final_Codon_Usage_Summary.csv'.")
    else:
        print("No data extracted.")