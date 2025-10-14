#!/bin/bash
set -e  # stop if any command fails

python scripts/clean_data.py -i data/processed/rat/full_data.csv -s rat -o data/processed/rat/cleaned_full_data.csv
python scripts/paper_blast.py -i data/processed/rat/cleaned_full_data.csv -o data/processed/rat/full_paper_blast.csv -s rat
python scripts/total_blast.py -i data/processed/rat/full_paper_blast.csv -s rat -o data/processed/rat/full_total_blast.csv
python scripts/align_to_full_seq.py -i data/processed/rat/full_total_blast.csv -s rat
