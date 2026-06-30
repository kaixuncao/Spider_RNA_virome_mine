import sys
import os
import math
import argparse
import csv
from Bio import SeqIO
from Bio.Data import CodonTable

# ==========================================
# 1. Core Algorithm: Weight Calculation and CAI Calculation
# ==========================================

def calculate_weights_from_counts(counts_dict):
    """
    Calculate the relative fitness weight (w).
    """
    genetic_code = CodonTable.unambiguous_dna_by_id[1]
    aa_to_codons = {}
    for codon, aa in genetic_code.forward_table.items():
        aa_to_codons.setdefault(aa, []).append(codon)
    for codon in genetic_code.stop_codons:
        aa_to_codons.setdefault('*', []).append(codon)

    weights = {}
    for aa, codon_list in aa_to_codons.items():
        current_counts = []
        valid_codons = []
        for c in codon_list:
            val = counts_dict.get(c)
            if val is None:
                val = counts_dict.get(c.replace('T', 'U'), 0)
            current_counts.append(val)
            valid_codons.append(c)
        if not current_counts: continue
        max_count = max(current_counts)
        for i, codon in enumerate(valid_codons):
            weights[codon] = current_counts[i] / max_count if max_count > 0 else 0.0
    return weights

def manual_cai_calculation(sequence, weights):
    """
    Calculate the CAI for a single sequence.
    Please note: if the original sequence is retained, it may contain a stop codon.
    Stop codons are not usually included in the weights dictionary; they are skipped at this stage and do not affect the CAI calculation.
    """
    seq = str(sequence).upper().replace('U', 'T')
    length = len(seq)
    if length % 3 != 0:
        seq = seq[:length - (length % 3)]
    if not seq: return None

    L = 0
    log_sum = 0.0
    
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        # If a stop codon or an unknown codon (not listed in the weight table) is encountered, simply skip it
        if codon in weights:
            w = weights[codon]
            if w > 0:
                log_sum += math.log(w)
                L += 1
    
    if L == 0: return None
    return math.exp(log_sum / L)

# ==========================================
# 2. Smart Extraction + Forced Length Retention
# ==========================================

def get_coding_length(seq):
    """The length of the sequence up to the first stop codon"""
    seq = str(seq).upper().replace('U', 'T')
    stop_codons = {'TAA', 'TAG', 'TGA'}
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if len(codon) < 3: return i 
        if codon in stop_codons: return i
    return len(seq) - (len(seq) % 3)

def extract_sequence_smart(full_seq, start, end):
    """
    Strategy:
    1. Calculate the length of the target coordinates (Target_Len).
    2. Run Frame Correction and ATG Optimization。
    3. Compare the final length with Target_Len.
       If the optimised sequence loses too much of its length (for example, < 50% of Target_Len),
       This indicates that too much information has been omitted in order to avoid the stop codon.
       At this point, a forced rollback to the original coordinate sequence (Frame 0, even if it contains a Stop) takes place.
    """
    s = start - 1
    e = end
    
    if s < 0 or e > len(full_seq): return "", "OutOfBounds"
    
    # Raw coordinate tiles
    raw_slice = str(full_seq[s:e]).upper()
    target_len = len(raw_slice) # Expected length
    
    if not raw_slice: return "", "Empty"

    # --- Phase A: Attempt to optimise the frame ---
    best_frame_seq = ""
    best_frame_len = -1
    best_offset = 0
    
    for offset in range(3):
        temp_seq = raw_slice[offset:]
        valid_len = get_coding_length(temp_seq)
        if valid_len > best_frame_len:
            best_frame_len = valid_len
            best_frame_seq = temp_seq[:valid_len]
            best_offset = offset
            
    base_seq = best_frame_seq
    current_status = f"Frame_{best_offset}"

    # --- Phase B: Attempting to optimise ATG ---
    final_seq = base_seq
    final_len = len(base_seq)
    
    search_idx = 0
    best_atg_seq = ""
    best_atg_len = -1
    
    while True:
        atg_idx = base_seq.find("ATG", search_idx)
        if atg_idx == -1: break
        if atg_idx % 3 == 0:
            candidate = base_seq[atg_idx:]
            c_len = len(candidate) 
            if c_len > best_atg_len:
                best_atg_len = c_len
                best_atg_seq = candidate
        search_idx = atg_idx + 1

    use_atg = False
    if best_atg_len > 0 and final_len > 0:
        ratio = best_atg_len / final_len
        # ATG must account for a certain proportion of the optimised sequence, and its absolute length must meet the required standard
        if (ratio >= 0.3 or final_len < 60) and best_atg_len > 6:
            use_atg = True

    if use_atg:
        final_seq = best_atg_seq
        current_status += "_ATG"
    
    # --- Phase C: [New] Length retention check ---
    # Calculate the length of the currently selected sequence
    current_len = len(final_seq)
    
    # Threshold setting: If the current length is less than 50% of the expected coordinate length
    # Explain that the optimisation process (due to the search for stop codons or ATG) discards too many sequences
    length_ratio = current_len / target_len if target_len > 0 else 0
    
    if length_ratio < 0.5:
        # Backtracking strategy: Return to the original sequence of coordinates (Frame 0), applying only simple pruning by multiples of 3
        # Even if there is a stop codon in the middle, it is still better than having only 84 bp left
        fallback_len = target_len - (target_len % 3)
        final_seq = raw_slice[:fallback_len]
        current_status = "Force_Raw_Length" # Marked as a mandatory rollback
    
    # The ‘Short’ flag is only set if the final sequence is genuinely short (and not due to a forced backtrack).
    elif len(final_seq) < 30:
        current_status += "_Short"

    return final_seq, current_status

# ==========================================
# 3. Data Loading
# ==========================================

def load_all_reference_tables(ref_dir):
    reference_data = {}
    if not os.path.exists(ref_dir):
        print(f"Error: Contents '{ref_dir}' Does not exist.")
        sys.exit(1)
    files = [f for f in os.listdir(ref_dir) if f.endswith(".txt")]
    if not files:
        print(f"Error: No .txt files were found in '{ref_dir}'.")
        sys.exit(1)
    print(f"Loading {len(files)} reference species lists...")
    for filename in files:
        filepath = os.path.join(ref_dir, filename)
        counts = {}
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("Codon"): continue
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        codon = parts[0].upper().replace('U', 'T')
                        try:
                            val = float(parts[1])
                            counts[codon] = val
                        except ValueError: pass
            weights = calculate_weights_from_counts(counts)
            clean_name = filename.replace(".gbff.txt", "").replace(".txt", "")
            reference_data[clean_name] = weights
        except Exception: pass
    print("The reference table has finished loading。")
    return reference_data

# ==========================================
# 4. Main process
# ==========================================

def main():
    parser = argparse.ArgumentParser(description="Batch CAI Calculator with Length Preservation")
    parser.add_argument("--fasta", required=True, help="Path to genomic FASTA.")
    parser.add_argument("--coords", required=True, help="Path to coordinates file.")
    parser.add_argument("--ref_dir", default="reference_tables", help="Directory of reference .txt files.")
    parser.add_argument("--output", default="CAI_Results.csv", help="Output CSV filename.")
    args = parser.parse_args()

    ref_weights_map = load_all_reference_tables(args.ref_dir)
    sorted_species_names = sorted(ref_weights_map.keys())
    
    print(f"Reading FASTA: {args.fasta} ...")
    if not os.path.exists(args.fasta):
        print("FASTA not found.")
        return
    seq_records = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))
    print(f"Loaded {len(seq_records)} sequences.")

    delimiter = ',' 
    has_header = False
    try:
        with open(args.coords, 'r') as f:
            line = f.readline()
            if not line: return
            if '\t' in line: delimiter = '\t'
            parts = line.strip().split(delimiter)
            try: int(parts[1])
            except: has_header = True
    except Exception as e:
        print(f"Error: {e}")
        return

    print(f"Processing: {args.coords} (Delimiter: {'Tab' if delimiter=='\\t' else 'Comma'})")

    match_count = 0
    miss_count = 0
    
    with open(args.coords, 'r') as infile, open(args.output, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter=delimiter)
        writer = csv.writer(outfile)
        writer.writerow(["Seq_ID", "Start", "End", "Status", "Length"] + sorted_species_names)

        if has_header: next(reader, None)

        for row in reader:
            if not row or len(row) < 3: continue
            row = [c.strip() for c in row]
            seq_id, start_str, end_str = row[0], row[1], row[2]
            try: start, end = int(start_str), int(end_str)
            except ValueError: continue

            if seq_id not in seq_records:
                miss_count += 1
                writer.writerow([seq_id, start, end, "ID_Not_Found", 0] + ["NA"]*len(sorted_species_names))
                continue

            full_seq = str(seq_records[seq_id].seq)
            target_seq, status = extract_sequence_smart(full_seq, start, end)

            if not target_seq or "Bound" in status:
                writer.writerow([seq_id, start, end, status, 0] + ["NA"]*len(sorted_species_names))
                continue

            cai_results = []
            for species in sorted_species_names:
                cai_val = manual_cai_calculation(target_seq, ref_weights_map[species])
                cai_results.append(f"{cai_val:.4f}" if cai_val is not None else "NA")

            writer.writerow([seq_id, start, end, status, len(target_seq)] + cai_results)
            match_count += 1
            if match_count % 1000 == 0: print(f"Processed {match_count} sequences...")

    print("-" * 50)
    print(f"Done! Processed {match_count} sequences.")
    print(f"ID Not Found: {miss_count}")
    print(f"Results saved to {args.output}")

if __name__ == "__main__":
    main()