# SLIDER (Scoring Sorting Screens with Laplace and Differential Rank)
from Lalit R. Patel

SLIDER identifies genes that are enriched or depleted in CRISPR screens that use fluorescence activated cell sorting (FACS).

Command Line:
```
  Rscript SLIDER_ScoreReplicates.R [Unsorted.cts] [Sorted.cts] [Name_Your_Replicate]
```

Input files (unsorted.cts, sorted.cts) have no header and four tab separated columns: sgRNA, additional_info, HUGO_GeneSymbol, counts. In the example below, the additional_info column indicates chromosomal position of sgRNA target site. Alternatively, this column can be used to store information about which guide pool the sgRNA comes from (for split pool libraries), whether the guides are targeting vs non-targeting, or anything else that might be useful. If you have no additional information for your sgRNA, fill column 2 with "none". Avoid using "NA", "na", or "NULL" in column 2.  For screens of a tiling library (i.e. enhancer mapping), assign bins every 5-6 sgRNA and list unique BinIDs or BinPositions in column 3 instead of gene symbol.

Example format of .cts file:
```
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
```

For each replicate scored using this script you will have six output files:
```
  Name_SLIDER_EnrichedGenes_ActiveGuides.csv
  Name_SLIDER_DepletedGenes_ActiveGuides.csv
  Name_SLIDER_NeutralGenes_NeutralGuides.csv 
  Name_SLIDER_sgRNA_scores_categorization.csv
  Name_SLIDER_GeneScores.csv
  Name_SLIDER_plots.pdf
```
	
If you have multiple replicates, first score each replicate using this command line script. Then open R or Rstudio, import the Name_SLIDER_GeneScores.csv output files generated for each of your replicates, and combine the GeneName and Zequiv columns from each output into a single dataframe using the merge function as follows:
```
  df_list <-list(Rep1_SLIDER_GeneScores[,c(2,8)],Rep2_SLIDER_GeneScores[,c(2,8)],…,RepN_SLIDER_GeneScores[,c(2,8)])
  dfmerged<-Reduce(function(x, y) merge(x, y, by=’GeneName’, all=TRUE), df_list)
```

Once you have a merged dataframe with Zequivs from your replicates, use colnames() to rename the columns of dfmerged to match the name of your replicates and store the number of replicates in the variable N:
```
  colnames(dfmerged)<-c(‘GeneName’,’Rep1’,’Rep2’,…’RepN’)
  N <- # of replicates
```
Use the following commands in R to calculate gene statistics across your replicates for a screen level score, pvalue, and FDR:
```
  dfmerged$Ns<- rowSums(!is.na(dfmerged[, c(2:(N+1))]))
  dfmerged$Zs = rowSums(dfmerged[,c(2:(N+1))],na.rm=TRUE) / sqrt(dfmerged$Ns)
  dfmerged$pvalue = 2*pnorm(q=abs(dfmerged$Zs), lower.tail=FALSE)
  dfmerged$FDR = p.adjust(dfmerged$pvalue, method = 'BH')
  Summary<-(dfmerged[order(-dfmerged$Zs),])
  View(Summary)
  write.csv(x=Summary, “ScreenSummary.csv”) 
```




