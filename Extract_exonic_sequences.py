#!/usr/bin/env python3



#This code use annotations from HG38.rna with exonic and intronic positions amd create a fasta file with the fasta sequences

######## Example : python extract_exonic_sequences.py PRKN_pos.rna hg38.fa PRKN_RNA.fasta


import sys
from Bio import SeqIO
from Bio.Seq import Seq
import re

def reverse_complement(seq):
    """Return reverse complement of DNA sequence"""
    return str(Seq(seq).reverse_complement())

def parse_annotations(annotation_string):
    """Parse the annotation string and return list of (element_name, start_pos) tuples"""
    elements = []
    parts = annotation_string.split(',')
    
    for part in parts:
        if ':' in part:
            element_name, pos = part.split(':')
            elements.append((element_name, int(pos)))
    
    return elements

def get_exonic_regions(elements):
    """Extract only exonic regions (starting with 'E') and calculate their coordinates"""
    exonic_regions = []
    
    for i, (element_name, start_pos) in enumerate(elements):
        if element_name.startswith('E'):
            # Find the end position (start of next element or use a default)
            if i + 1 < len(elements):
                end_pos = elements[i + 1][1]
            else:
                # For the last element, we need to estimate the end
                # This is a limitation - you might want to provide end coordinates
                end_pos = start_pos + 20  # Default assumption
                print(f"Warning: No end coordinate for last element {element_name}, using +20bp")
            
            exonic_regions.append((element_name, start_pos, end_pos))
    
    return exonic_regions

def extract_sequence_from_fasta(fasta_file, chromosome, start, end):
    """Extract sequence from FASTA file given coordinates (1-based)"""
    try:
        # Load the chromosome sequence
        if chromosome not in fasta_file:
            raise KeyError(f"Chromosome {chromosome} not found in FASTA file")
        
        chr_seq = fasta_file[chromosome]
        # Convert to 0-based indexing for Python
        sequence = str(chr_seq.seq[start-1:end])
        return sequence
    except Exception as e:
        print(f"Error extracting sequence {chromosome}:{start}-{end}: {e}")
        return ""

def process_rna_file(input_file, hg38_fasta_path, output_file):
    """Main function to process the RNA annotation file"""
    
    print("Loading hg38 FASTA file...")
    # Load hg38 FASTA file (this might take a while)
    fasta_sequences = SeqIO.to_dict(SeqIO.parse(hg38_fasta_path, "fasta"))
    print(f"Loaded {len(fasta_sequences)} sequences from hg38")
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line_num, line in enumerate(infile, 1):
            line = line.strip()
            if not line:
                continue
                
            try:
                # Parse the tab-separated line
                parts = line.split('\t')
                if len(parts) < 6:
                    print(f"Warning: Line {line_num} has fewer than 6 columns, skipping")
                    continue
                
                rna_id = parts[0]
                chromosome = parts[1]
                strand = parts[4]
                annotations = parts[5]
                
                print(f"Processing {rna_id} on {chromosome} ({strand})")
                
                # Parse annotations
                elements = parse_annotations(annotations)
                
                # Get exonic regions
                exonic_regions = get_exonic_regions(elements)
                
                if not exonic_regions:
                    print(f"Warning: No exonic regions found for {rna_id}")
                    continue
                
                # Extract sequences for each exonic region
                concatenated_sequence = ""
                
                for element_name, start_pos, end_pos in exonic_regions:
                    sequence = extract_sequence_from_fasta(fasta_sequences, chromosome, start_pos, end_pos)
                    if sequence:
                        concatenated_sequence += sequence
                        print(f"  {element_name}: {start_pos}-{end_pos} ({len(sequence)}bp)")
                
                # Handle strand orientation
                if strand == '-':
                    concatenated_sequence = reverse_complement(concatenated_sequence)
                    print(f"  Reverse complemented for minus strand")
                
                # Write to output FASTA
                if concatenated_sequence:
                    outfile.write(f">{rna_id}\n")
                    # Write sequence in 80-character lines
                    for i in range(0, len(concatenated_sequence), 80):
                        outfile.write(concatenated_sequence[i:i+80] + '\n')
                    
                    print(f"  Final sequence length: {len(concatenated_sequence)}bp")
                else:
                    print(f"Warning: No sequence extracted for {rna_id}")
                    
            except Exception as e:
                print(f"Error processing line {line_num}: {e}")
                continue
    
    print(f"Done! Output written to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extract_exonic_sequences.py <input_file> <hg38_fasta> <output_fasta>")
        print("Example: python extract_exonic_sequences.py annotations.txt hg38.fa output.fasta")
        sys.exit(1)
    
    input_file = sys.argv[1]
    hg38_fasta_path = sys.argv[2]
    output_file = sys.argv[3]
    
    process_rna_file(input_file, hg38_fasta_path, output_file)
