# Rat Data Decisions/Cleaning/Analysis

## Sep 14 2025
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