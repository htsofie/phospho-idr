Scripts used for data cleaning, analysing and workflow

### Workflow:
gen_data.py -> clean_data.py -> map_to_uniprot.py -> get_full_seq.py -> align_to_full_seq.py

*need to write down CLI commands to use scripts and their pathway.
* this will probably be completed with nextflow at somepoint

gen_data.py: 
    Generates either full or sample data set in either csv or parquet file. Need to specify which config you want.
    This script is run through the command line.

align_to_full_seq.py:

'''
#command line prompts to run script
    python scripts/align_to_full_seq.py --input data/processed/mouse/cleaned_sample_mapped_sequences.csv --species mouse

    python scripts/align_to_full_seq.py --input data/processed/rat/cleaned_sample_mapped_sequences.csv --species rat
'''