# ApoScraper

Uses pypdb to scrape the 'related proteins' from the PDB for a specified PDB code and returns PDB codes deemed to be apo structures

Note: Needs refinement for how related proteins are chosen. Needs to be changed to a specified sequence similarity threshold.

Dependencies

1. pypdb
    ``` pip install pypdb```
 
Running the apo_scraper.py script

``` python apo_scraper {PDB_ID} ```

- PDB_ID = any 4 letter PDB_ID code
