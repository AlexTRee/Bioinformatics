Regeneron Legacy Projects Visualization
==========================================

## Description

- VelociGene, that uses targeting vectors based on bacterial artificial chromosomes (BACs), can precisely replace the gene of interest with a reporter

- Since 2012, VelociWeb parses design sequence elements (recombination oligos, deletions and TaqMan assays) to the Regn-UCSC genome browsers via a weekly cron job that extracts information from the Targeting Genbank files from every new MAID

    The 1st step is to extract desired information from the Jweb database
    The 2nd step is to find restriction enzyme cutting sites and remove primers’ tails before BLAT search
    The 3rd step is to use customized BLAT search to determine the sequence coordinates
      - BLAT server has been set up with parameters that can find shorter sequences (min=15bps, less than default 20bps)
    The 4th step is to filter BLAT results (HUR-HUF < 600, TUR-TUF < 500, No overlap)
    The 5th step is to calculate coordinates for deletions 
      - The Deletion is defined by the region between HUR and HDF
    The last step is to convert data into .BED file format for UCSC Genome Browser visualization
      - chr7	16822253	16822272	100HDF	100	-

## Final Result:
# 98.3% of legacy MAIDs are successfully matched with deletion.

