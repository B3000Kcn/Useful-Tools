#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: abricate2gbk.py
Description: 
    A specialized tool to convert Abricate mass screening results (TSV) 
    into standard GenBank format (.gbk).
    
    This script addresses the limitation of raw TSV output by mapping 
    resistance gene annotations to BioPython SeqFeature objects, enabling 
    direct visualization in software like SnapGene and Artemis.

    Key Features:
    - Robust CSV/TSV parsing.
    - Auto-calculation of feature locations (1-based to 0-based conversion).
    - Metadata injection (Gene, Product, %Identity, %Coverage).
    - Change the type in line 123 to alter the type of annotation (e.g., Gene, tRNA, etc.)

Dependencies:
    - biopython>=1.81
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
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

# --- Main Logic ---

if len(sys.argv) != 4:
    sys.stderr.write(f"Usage: {sys.argv[0]} input.fasta abricate.tsv output.gbk\n")
    sys.exit(1)

fasta_file = sys.argv[1]
abricate_file = sys.argv[2]
out_file = sys.argv[3]

# 1. Read FASTA -> dict: seq_id -> SeqRecord
seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

# 2. Prepare a list of features for each SeqRecord
features_by_id = defaultdict(list)

# 3. Read Abricate TSV
with open(abricate_file, "r", newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    
    # Ensure required columns exist
    required = ["SEQUENCE", "START", "END", "STRAND", "GENE", "DATABASE"]
    for r in required:
        if r not in reader.fieldnames:
            sys.stderr.write(f"[ERROR] Missing column in abricate file: {r}\n")
            sys.exit(1)

    # Check for optional columns
    has_cov = "%COVERAGE" in reader.fieldnames
    has_ident = "%IDENTITY" in reader.fieldnames
    has_product = "PRODUCT" in reader.fieldnames
    has_res = "RESISTANCE" in reader.fieldnames

    count = 0
    for row in reader:
        seq_id = row["SEQUENCE"]
        if seq_id not in seq_dict:
            # Skip if sequence ID in Abricate is not found in FASTA
            continue

        try:
            start = int(row["START"])
            end = int(row["END"])
        except ValueError:
            # Malformed line (invalid start/end), skip
            continue

        strand_char = row["STRAND"]
        strand = 1
        if strand_char == "-":
            strand = -1

        gene = row["GENE"]
        db = row["DATABASE"]

        # Handle optional fields: product, coverage, identity, resistance
        product = row["PRODUCT"] if has_product and row["PRODUCT"] else gene
        cov = row["%COVERAGE"] if has_cov and row["%COVERAGE"] else ""
        ident = row["%IDENTITY"] if has_ident and row["%IDENTITY"] else ""
        resistance = row["RESISTANCE"] if has_res and row["RESISTANCE"] else ""

        # Biopython FeatureLocation is 0-based, [start, end)
        location = FeatureLocation(start - 1, end, strand=strand)

        qualifiers = {
            "gene": [gene],
            "product": [product],
        }

        notes = []
        if cov:
            notes.append(f"Coverage={cov}%")
        if ident:
            notes.append(f"Identity={ident}%")
        if db:
            notes.append(f"Database={db}")
        if resistance:
            notes.append(f"Resistance={resistance}")

        if notes:
            qualifiers["note"] = ["; ".join(notes)]

        # Create the feature (Type is set to 'CDS')
        feat = SeqFeature(location=location, type="CDS", qualifiers=qualifiers)
        features_by_id[seq_id].append(feat)
        count += 1

sys.stderr.write(f"[INFO] Loaded {count} abricate hits.\n")

# 4. Attach features to SeqRecord and write to GenBank
records_out = []
for seq_id, rec in seq_dict.items():
    # Create a new SeqRecord to avoid modifying the original one
    new_rec = SeqRecord(
        rec.seq,
        id=rec.id,
        name=rec.name,
        description=rec.description
    )
    # GenBank requires basic annotations to avoid warnings in some software (e.g., SnapGene)
    new_rec.annotations["molecule_type"] = "DNA"
    new_rec.annotations["topology"] = "circular"  # Assumed circular (standard for bacterial plasmids)
    new_rec.annotations["data_file_division"] = "BCT"

    new_rec.features = features_by_id.get(seq_id, [])
    records_out.append(new_rec)

with open(out_file, "w") as out_handle:
    SeqIO.write(records_out, out_handle, "genbank")

sys.stderr.write(f"[INFO] Wrote {len(records_out)} records to {out_file}\n")