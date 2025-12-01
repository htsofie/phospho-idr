#!/bin/bash
set -e  # stop if any command fails

python scripts/clean_data.py -i data/processed/mouse/full_data.csv -s mouse -o data/processed/mouse/cleaned_full_data.csv
python scripts/paper_blast.py -i data/processed/mouse/cleaned_full_data.csv -o data/processed/mouse/full_paper_blast.csv -s mouse
python scripts/total_blast.py -i data/processed/mouse/full_paper_blast.csv -s mouse -o data/processed/mouse/full_total_blast.csv
python scripts/align_to_full_seq.py -i data/processed/mouse/full_total_blast.csv -s mouse
