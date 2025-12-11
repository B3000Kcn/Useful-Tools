#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: abricate2gff.py
Description: 
    A specialized tool to convert Abricate mass screening results (TSV) 
    into standard GFF3 format (.gff).
    
    This is the sibling tool to abricate2gbk, designed for pipelines 
    that require GFF3 input (e.g., Roary, JBrowse).

    Key Features:
    - 1-based coordinate system preservation (No conversion needed).
    - Generates standard GFF3 headers including ##sequence-region.
    - Encodes Abricate metadata (Identity, Coverage, Database) into GFF3 attributes.
    - Sanitizes attribute strings to ensure GFF3 compatibility.

    How to use: abricate2gff.py input.fasta abricate.tsv output.gff

Dependencies:
    - biopython>=1.81 (for parsing FASTA lengths)
"""

# --- Metadata ---
__author__      = "B3000Kcn"
__credits__     = ["DBL1F7E5"]
__copyright__   = "Copyright 2025, B3000Kcn"
__license__     = "MIT"
__version__     = "1.0.0"

# --- Imports ---
import sys
import csv
import urllib.parse
from Bio import SeqIO

# --- Main Logic ---

def sanitize_gff_attr(key, value):
    """
    Sanitize values for GFF3 attributes.
    Quotes are percent-encoded if necessary, semicolons/equals are escaped.
    """
    if not value:
        return ""
    # Simple replacement for common forbidden chars in GFF3 attributes
    # Or strict percent encoding: urllib.parse.quote(str(value))
    # Here we keep it readable but safe-ish
    val_str = str(value).replace(";", "%3B").replace("=", "%3D").replace("&", "%26")
    return val_str

def main():
    if len(sys.argv) != 4:
        sys.stderr.write(f"Usage: {sys.argv[0]} input.fasta abricate.tsv output.gff\n")
        sys.exit(1)

    fasta_file = sys.argv[1]
    abricate_file = sys.argv[2]
    out_file = sys.argv[3]

    # 1. Read FASTA to get sequence lengths for ##sequence-region
    #    We don't need to load the whole seq into memory, just ID and length.
    seq_lengths = {}
    try:
        # parsing with "fasta" is generator-based, efficient
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq_lengths[record.id] = len(record.seq)
        sys.stderr.write(f"[INFO] Read {len(seq_lengths)} sequences from FASTA.\n")
    except Exception as e:
        sys.stderr.write(f"[ERROR] Failed to read FASTA: {e}\n")
        sys.exit(1)

    # 2. Process Abricate TSV and buffer features
    gff_lines = []
    
    with open(abricate_file, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        
        # Validation
        required = ["SEQUENCE", "START", "END", "STRAND", "GENE", "DATABASE"]
        for r in required:
            if r not in reader.fieldnames:
                sys.stderr.write(f"[ERROR] Missing column: {r}\n")
                sys.exit(1)

        has_cov = "%COVERAGE" in reader.fieldnames
        has_ident = "%IDENTITY" in reader.fieldnames
        has_product = "PRODUCT" in reader.fieldnames
        has_res = "RESISTANCE" in reader.fieldnames
        has_accession = "ACCESSION" in reader.fieldnames

        count = 0
        for i, row in enumerate(reader):
            seq_id = row["SEQUENCE"]
            
            # Skip if sequence not in FASTA (optional check, but good for consistency)
            if seq_id not in seq_lengths:
                continue

            try:
                start = int(row["START"])
                end = int(row["END"])
            except ValueError:
                continue

            # GFF3 Strand: '+', '-', or '.'
            strand = row["STRAND"]
            if strand not in ['+', '-']:
                strand = '.'

            # Source & Type
            source = row["DATABASE"] if row["DATABASE"] else "Abricate"
            feature_type = "gene" # or CDS, region, match... 'gene' is safest for viewing

            # Score: Using %IDENTITY as score is common practice, or '.'
            score = row["%IDENTITY"] if has_ident else "."

            # Attributes
            # ID is mandatory for good GFF3
            feat_id = f"abricate_{i+1}"
            gene = row["GENE"]
            
            attrs = [f"ID={feat_id}"]
            attrs.append(f"Name={sanitize_gff_attr('Name', gene)}")
            
            if has_product and row["PRODUCT"]:
                attrs.append(f"product={sanitize_gff_attr('product', row['PRODUCT'])}")
            
            if has_ident and row["%IDENTITY"]:
                attrs.append(f"identity={row['%IDENTITY']}") # Numbers don't need sanitizing
            
            if has_cov and row["%COVERAGE"]:
                attrs.append(f"coverage={row['%COVERAGE']}")
            
            if has_res and row["RESISTANCE"]:
                 attrs.append(f"resistance={sanitize_gff_attr('resistance', row['RESISTANCE'])}")

            if has_accession and row["ACCESSION"]:
                 attrs.append(f"accession={sanitize_gff_attr('accession', row['ACCESSION'])}")

            attributes_str = ";".join(attrs)

            # Construct GFF line (9 columns)
            # 1. seqid
            # 2. source
            # 3. type
            # 4. start
            # 5. end
            # 6. score
            # 7. strand
            # 8. phase (CDS need 0/1/2, genes use '.')
            # 9. attributes
            line = f"{seq_id}\t{source}\t{feature_type}\t{start}\t{end}\t{score}\t{strand}\t.\t{attributes_str}"
            gff_lines.append(line)
            count += 1

    sys.stderr.write(f"[INFO] Processed {count} annotations.\n")

    # 3. Write GFF3 File
    with open(out_file, "w") as out:
        # Header
        out.write("##gff-version 3\n")
        
        # Sequence-region directives (Standard GFF3 practice)
        for seq_id, length in seq_lengths.items():
            out.write(f"##sequence-region {seq_id} 1 {length}\n")
        
        # Features
        for line in gff_lines:
            out.write(line + "\n")
        
        # Optional: FASTA at the end (some tools require this)
        # out.write("##FASTA\n") 
        # But we skip it here to keep file size small; it's separate usually.

    sys.stderr.write(f"[INFO] Done. Output written to {out_file}\n")

if __name__ == "__main__":
    main()
