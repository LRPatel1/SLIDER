# SLIDER
Scoring Soriting Screens with Laplace and Differential Rank

SLIDER identifies genes that are enriched or depleted in CRISPR screens using fluorescence activated cell sorting (FACS).


Command Line:
  Rscript SLIDER_ScoreReplicates.R [Unsorted.cts] [Sorted.cts] [Name_Your_Replicate]
  
Input files (unsorted.cts, sorted.cts) have no header and four columns tab separated columns: sgRNA, additional_info, GeneName, counts. In the example below, the additional_info column indicates chromosomal position of sgRNA target site. You can also use additional info to indicate which guide pool (for split pool libraries), whether the guides are targeting vs non-targeting, or anything else. If you have no additional information for your sgRNA, fill column 2 with ‘NA’. 

Example format of .cts file:
GACAAAGCGCGCCAGCGCGG  chr1:220960306-220960325-exon01 43525 267
GAAGTCACAGCCCTACCGCC  chr1:220970054-220970073+exon03 43525 2487
CATGAGCAGGAAGGAACCGC  chr1:220978581-220978600+exon06 43525 92
ATGGATACTGTACCTTCCGG  chr4:164506888-164506907+exon06 43525 3291
AGGCAAGCAGCCCAACAACA  chr4:164534528-164534547-exon05 43525 1881
GATCGACTTGCAGATCGCCC  chr4:164621981-164622000+exon04 43525 607
GCGAACACCCGAGATCTCAG  chr4:164775208-164775227-exon03 43525 222
ATGGAGATGAGCACGAGGCG  chr1:220928344-220928363-exon02 43526 41
GAGCGGGCAGTAGTCTGGGT  chr1:220936259-220936278-exon04 43526 718
GGTATTCAAATCTACCAGGG  chr1:220936298-220936317-exon04 43526 32



