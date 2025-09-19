# General Decisions - data cleaning/workflow/pipeline

## Rat Localization Filtering

### Sep 14 2025
Which collumns to keep for data standardization.
Looking at localization probability. It is the probability that the phosphorylation site identified is highly localized. In notebooks/search_data.py I looked at the max and min values of localization prob across the dataset to determine if I had to remove any low localization values. 
    Max - 1.00
    Min - 0.745606
  According to the research paper the data is retrieved from, they only kept localized scores >0.75 however it seems that there are some scores less - 0.745. In notebooks I will try to identify how many are less than the hard cutoff and remove those. 
    Number of rows with localization prob < 0.75 = 11
So I need to remove these 11 rows from the dataset
First though, I want to check what their values are.
    1741     0.748312
    1885     0.749659
    2962     0.749088
    3943     0.745606
    7013     0.746718
    8324     0.747437
    8396     0.746815
    9252     0.749825
    14414    0.749021
    20591    0.749738
    21383    0.748941   
To remove these later will add a filter to config for localizaiton prob < 0.75 

## Protein ID mapping to Uniprot to get full sequence

### Sept 15 2025
    Trying to decide how to map protein ids from datasets to proteins in uniprot and how to ensure that the phosphosite and sequence window from the dataset match to the full sequence with the correct position number.

    Debated looking back timewise to get a protein sequence upload that is consistent around the same time of the data set generation. However, there is no way to really tell that the full sequence uploaded is from that data set. I think this is also a bit extra and a fragile method. Instead I could just try to map the sequence window to the canonical sequence or most recent sequence in uniprot for the isoform or protein and then see if it matches somewhere and adjust position value for new sequence and pull full sequence for the protein. 

    It would look something like this:
    swiss-prot/uniprot/ensembl id -> uniprot id -> most recent full sequence -> search for area that matches sequence window -> if multiple return the one closest to the position of phosphorylation documented in the datasets -> if none match return any that have only 1-4 AAs different -> if still none match then return which lines don't match or save in new data set. -> extract full sequence from best match -> add full sequence to data and edit position of phosphorylation in sequence to be the correct one. 

### Sept 16th 2025
    After mapping IDs to uniprot with to confirm/retrieve corresponding uniprot id to make sure it still exists, some proteins have been removed due to redundancy/have been invalidated, I extract full sequence with that uniprot id for either isoform - if specified - or canonical sequence. 

    To confirm this is done correctly I ran the script on my sample data csvs and am manually confirming in the uniprot database that the id matches the correct protein. I am cross referencing this with the protein name/description in the raw data set. 

    Then I will align the phosphorylation sequence window with the full protein sequence.

### Sept 17th 2025
    Looking at alignment of sequence window to full protein, I see a couple proteins where the alignment is very far off from the expected position of the phosphorylation. Some of them, the proteins identified, even have an expected position that is larger than the uniprot id mapped protein itself. I compared this to the IDs of the proteins in the raw data set and it is consistent. In uniprot showing a protein length of 123 AA or 198 but the phosphorylation position is 1293 or something which can't be possible.

    In get_full_seq.py I am flagging any rows where the sequence length number is smaller than the position number of the phosphorylation. These rows, depending on how many proteins there are, I can redirect to a different ID and sequence to align to. 