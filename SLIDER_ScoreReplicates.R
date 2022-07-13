#!/usr/bin/env SLIDER_ScoreReplicates.R
#!/usr/local/bin/ Rscript

args = commandArgs(trailingOnly = TRUE)
if (length(args)<3) {
  stop('Missing one or more input arguments. Provide (counts file for unsorted sample), (counts file for sorted sample), and (a name you will use to identify this analysis)', call.=FALSE)
}
if (length(args)>3) {
  stop('You have too many input arguments. Provide (counts file for unsorted sample), (counts file for sorted sample), and (a name you will use to identify this analysis)', call.=FALSE)
}

#Read in supplied counts data and title/label for this sort
Unsort <- read.table(args[1], header=FALSE, row.names=1, quote='\'', comment.char='', stringsAsFactors = FALSE)
Sort <- read.table(args[2], header=FALSE, row.names=1, quote='\'', comment.char='', stringsAsFactors = FALSE)
StudyName <- as.character(args[3])

#Create a dataframe for read counts at each guide present in Sorted and Unsorted samples
Data <- merge(Unsort, Sort, by=0, all.x=TRUE)
Data$V2.y<-NULL
Data$V3.y<-NULL

#Assigning a value of 0 counts to sgRNA present as NAs (these are guides present in the Unsorted or Sorted Sample but not captured by NGS of the other)
Data[is.na(Data)]<-0

#Assigning Labels to Values 
Data$AdditionalGuideInfo<-Data$V2.x
Data$V2.x<-NULL
Data$GeneName<-Data$V3.x
Data$V3.x<-NULL
Data$UnsortCounts<-Data$V4.x
Data$SortCounts<-Data$V4.y
Data$V4.x<-NULL
Data$V4.y<-NULL

#Normalizing Counts for Read Depth Variation (Counts Per Million) and log transforming using log2(CPM+2) pseudocount for plotting
UnsortSizeFactor<-(sum(Data$UnsortCounts))/1000000
SortSizeFactor<-(sum(Data$SortCounts))/1000000
Data$UnsortNormCounts<-(Data$UnsortCounts/UnsortSizeFactor)
Data$SortNormCounts<-(Data$SortCounts/SortSizeFactor)
Data$Unsort_log2PC <- log2(Data$UnsortNormCounts+2)
Data$Sort_log2PC <- log2(Data$SortNormCounts+2)

#Assigning sgRNA Ranks for Each Sample Based on Raw Read Counts and Calculating Difference in Rank for Statistical Testing
Data$Unsort_Rank<-rank(x=Data$UnsortCounts)
Data$Sort_Rank<-rank(x=Data$SortCounts)
Data$Unsort_NormRank<-Data$Unsort_Rank/length(Data$Unsort_Rank)
Data$Sort_NormRank<-Data$Sort_Rank/length(Data$Sort_Rank)
Data$SortVsUnsort_NormRankDiff<-Data$Sort_NormRank-Data$Unsort_NormRank

#Dividing Data Set into sgRNAs that increased in Rank (positive rank diff) and sgRNAs that decreased or stayed the same (rank shift <=0)
DataDown<-Data[which(Data$SortVsUnsort_NormRankDiff <= 0),]
DataUp<-Data[which(Data$SortVsUnsort_NormRankDiff > 0),]

#Computing Local Estimates for the Diversity Parameter (b) for a Laplace Distribution for the Normalized Change in sgRNA Ranks (Sort vs Unsort) by  
#sampling 1000 sgRNA closest in rank to the current sgRNA within the unsorted dataset and defining an expected value of 0 for a null change in rank.
#Using the maximum likelihood estimator to calculate b. Performing estimate for two half-functions, one for guides that enriched after sorting (rank  
#shift>0) and one for guides that did not (rank shift<=0). Applying distribution with locally estimated directional diversity to compute sgRNA pvalue.
Laplace_pval <- function(a,mu,b) {
  p<-0.5*(exp(-1*abs(a-mu)/b))
  p_twotail<-2*p
  return (p_twotail)
}
center = 0
diversity = sum(abs(DataUp$SortVsUnsort_NormRankDiff))/length(DataUp$SortVsUnsort_NormRankDiff)
DataUp_Ordered<-DataUp[order(DataUp$UnsortNormCounts),]
DataUp_Ordered$idx<-seq(from=1,to=(length(DataUp$SortVsUnsort_NormRankDiff)))
DataUp_Ordered$minRange<-DataUp_Ordered$idx - 499
DataUp_Ordered$maxRange <-DataUp_Ordered$idx + 500
DataUp_Ordered$minRange[which(DataUp_Ordered$idx < 500)] <- 1
DataUp_Ordered$maxRange[which(DataUp_Ordered$idx < 500)] <- 1000
DataUp_Ordered$minRange[which(DataUp_Ordered$idx > length(DataUp_Ordered$SortVsUnsort_NormRankDiff)-500)]<-length(DataUp_Ordered$SortVsUnsort_NormRankDiff)-999
DataUp_Ordered$maxRange[which(DataUp_Ordered$idx > length(DataUp_Ordered$SortVsUnsort_NormRankDiff)-500)]<-length(DataUp_Ordered$SortVsUnsort_NormRankDiff)
i=0
while(i<length(DataUp_Ordered$idx)){
  i=i+1
  Nearest1K<-DataUp_Ordered$SortVsUnsort_NormRankDiff[which(DataUp_Ordered$idx>=DataUp_Ordered$minRange[i] & DataUp_Ordered$idx<=DataUp_Ordered$maxRange[i])]
  DataUp_Ordered$b[i]<-sum(abs(Nearest1K - 0))/length(Nearest1K)
}
DataUp_Ordered$sgRNA_pvalue<- Laplace_pval(DataUp_Ordered$SortVsUnsort_NormRankDiff,mu=0,DataUp_Ordered$b)
DataUp_Ordered$idx<-NULL
DataUp_Ordered$minRange<-NULL
DataUp_Ordered$maxRange<-NULL

center = 0
diversity = sum(abs(DataDown$SortVsUnsort_NormRankDiff))/length(DataDown$SortVsUnsort_NormRankDiff)
DataDown_Ordered<-DataDown[order(DataDown$UnsortNormCounts),]
DataDown_Ordered$idx<-seq(from=1,to=(length(DataDown$SortVsUnsort_NormRankDiff)))
DataDown_Ordered$minRange<-DataDown_Ordered$idx - 499
DataDown_Ordered$maxRange <-DataDown_Ordered$idx + 500
DataDown_Ordered$minRange[which(DataDown_Ordered$idx < 500)] <- 1
DataDown_Ordered$maxRange[which(DataDown_Ordered$idx < 500)] <- 1000
DataDown_Ordered$minRange[which(DataDown_Ordered$idx > length(DataDown_Ordered$SortVsUnsort_NormRankDiff)-500)]<-length(DataDown_Ordered$SortVsUnsort_NormRankDiff)-999
DataDown_Ordered$maxRange[which(DataDown_Ordered$idx > length(DataDown_Ordered$SortVsUnsort_NormRankDiff)-500)]<-length(DataDown_Ordered$SortVsUnsort_NormRankDiff)
i=0
while(i<length(DataDown_Ordered$idx)){
  i=i+1
  Nearest1K<-DataDown_Ordered$SortVsUnsort_NormRankDiff[which(DataDown_Ordered$idx>=DataDown_Ordered$minRange[i] & DataDown_Ordered$idx<=DataDown_Ordered$maxRange[i])]
  DataDown_Ordered$b[i]<-sum(abs(Nearest1K - 0))/length(Nearest1K)
}
DataDown_Ordered$sgRNA_pvalue<- Laplace_pval(DataDown_Ordered$SortVsUnsort_NormRankDiff,mu=0,DataDown_Ordered$b)
DataDown_Ordered$idx<-NULL
DataDown_Ordered$minRange<-NULL
DataDown_Ordered$maxRange<-NULL
Data_Ordered<-rbind(DataUp_Ordered,DataDown_Ordered)

#An sgRNA that had no change in rank will recieve a p-value of 1 from the laplace function above. This will lead to logits that are infinite in steps below. Replacing
#these p-values with 0.99999995 to prevent 0 values from impairing calculations later in this script. 
Data_Ordered$sgRNA_pvalue[which(Data_Ordered$SortVsUnsort_NormRankDiff==0)]<-0.99999995

#The next step is a filter. Not all sgRNA efficiently knockout/activate/repress their targets in every cell line. Some sgRNAs also have cell-line specific off-target
#activity. This causes a minority of sgRNA designed to target a gene to behave differently from sgRNA that successfully target the gene. In a screen where evidence is 
#combined across independant sgRNA to nominate genes, retaining these sgRNA reduces sensitivity. To overcomes this, a filter is implemented to exclude sgRNA that are  
#discordant in their behavior from that of the predominate behavior of sgRNAs designed against a specific gene.
#
#Filtering Steps: 
##Identify Genes Where Fisher Test Finds the Fraction of sgRNAs that Increase in Rank at a Gene is Greater than the Fraction Across the Library. Label as UpGene.
#Identify Genes Where Fisher Test Finds the Fraction of sgRNAs that Decrease in Rank at a Gene is Greater than the Fraction Across the Library. Label as DownGene.
#Identify Genes Where Fisher Tests Found No Such Difference for Either Fraction. Label as Neutral. Genes categorized by these tests as both UpGenes and DownGenes  
#are relabeled as Neutral on the grounds that sorting lead to the sgRNAs for that gene becoming dispersed but without any prefferance for direction. Using a 
#permissive Fisher's Exact pvalue of 0.1 for this step. Feel free to tune this filtering parameter for the needs of your library, screen, and dataset. Once genes
#are labeled, those sgRNAs that are discordant with the genes category label are excluded from further analysis.

#Categorizing sgRNA with p<0.2 as active and enriched (positive rank change) or depleted (negative rank change). sgRNA unable to acheive p<0.2 after sorting are
#are categorized as neutral.
Data_Ordered$sgRNA_Up <-0
Data_Ordered$sgRNA_Up[which(Data_Ordered$sgRNA_pvalue<0.2 & Data_Ordered$SortVsUnsort_NormRankDiff>0)]<-1
Data_Ordered$sgRNA_Down <-0
Data_Ordered$sgRNA_Down[which(Data_Ordered$sgRNA_pvalue<0.2 & Data_Ordered$SortVsUnsort_NormRankDiff<0)]<-1
Data_Ordered$sgRNA_Neutral <-1
Data_Ordered$sgRNA_Neutral[which(Data_Ordered$sgRNA_Up==1)]<-0
Data_Ordered$sgRNA_Neutral[which(Data_Ordered$sgRNA_Down==1)]<-0

#Performing Fisher's Exact Test Comparing the Fraction of sgRNAs for Gene-i that Increase in Rank to the Fraction of sgRNAs that Do So Across the Whole Library. 
Gene_sgRNAcount<-aggregate(Data_Ordered$sgRNA_Up, by=list(Data_Ordered$GeneName),FUN=length)
Gene_sgRNAup<-aggregate(Data_Ordered$sgRNA_Up, by=list(Data_Ordered$GeneName),FUN=sum)
UpGeneSummary<-merge(Gene_sgRNAcount, Gene_sgRNAup, by='Group.1')
UpGeneSummary$GeneName<-UpGeneSummary$Group.1
UpGeneSummary$NsgRNA<-UpGeneSummary$x.x
UpGeneSummary$NsgRNAup<-UpGeneSummary$x.y
LibraryUp<-sum(Data_Ordered$sgRNA_Up)
LibraryTotal<-length(Data_Ordered$sgRNA_Up)
FisherExactPvalue<-Vectorize(function(NumUp,NumTotal,LibUp,LibTotal) {
  out<-fisher.test(matrix(c(NumUp, NumTotal-NumUp, LibUp, LibTotal-LibUp), ncol=2))
  return(out$p.value)
})
UpGeneSummary$EnrichedPvalue<-FisherExactPvalue(UpGeneSummary$NsgRNAup,UpGeneSummary$NsgRNA,LibraryUp,LibraryTotal)
UpGeneSummary$Group.1<-NULL
UpGeneSummary$x.x<-NULL
UpGeneSummary$x.y<-NULL

#Performing Fisher's Exact Test Comparing the Fraction of sgRNAs for Gene-i that Decrease in Rank to the Fraction of sgRNAs that Do So Across the Whole Library. 
Gene_sgRNAdown<-aggregate(Data_Ordered$sgRNA_Down, by=list(Data_Ordered$GeneName),FUN=sum)
DownGeneSummary<-merge(Gene_sgRNAcount, Gene_sgRNAdown, by='Group.1')
DownGeneSummary$GeneName<-DownGeneSummary$Group.1
DownGeneSummary$NsgRNA<-DownGeneSummary$x.x
DownGeneSummary$NsgRNAdown<-DownGeneSummary$x.y
LibraryDown<-sum(Data_Ordered$sgRNA_Down)
LibraryTotal<-length(Data_Ordered$sgRNA_Up)
DownGeneSummary$DepletedPvalue<-FisherExactPvalue(DownGeneSummary$NsgRNAdown,DownGeneSummary$NsgRNA,LibraryDown,LibraryTotal)
DownGeneSummary$Group.1<-NULL
DownGeneSummary$x.x<-NULL
DownGeneSummary$x.y<-NULL

UpGenes<-UpGeneSummary[which(UpGeneSummary$EnrichedPvalue<0.1),]
DownGenes<-DownGeneSummary[which(DownGeneSummary$DepletedPvalue<0.1),]
BothSets<-UpGenes[which(UpGenes$GeneName %in% DownGenes$GeneName),]
Enriched<-UpGenes[which(!(UpGenes$GeneName %in% BothSets$GeneName)),]
Depleted<-DownGenes[which(!(DownGenes$GeneName %in% BothSets$GeneName)),]

Gene_sgRNAneut<-aggregate(Data_Ordered$sgRNA_Neutral, by=list(Data_Ordered$GeneName),FUN=sum)
NeutGeneSummary<-merge(Gene_sgRNAcount, Gene_sgRNAneut, by='Group.1')
NeutGeneSummary$GeneName<-NeutGeneSummary$Group.1
NeutGeneSummary$NsgRNA<-NeutGeneSummary$x.x
NeutGeneSummary$NsgRNAneutral<-NeutGeneSummary$x.y
NeutGeneSummary$Group.1<-NULL
NeutGeneSummary$x.x<-NULL
NeutGeneSummary$x.y<-NULL
Neutral<-NeutGeneSummary[which(!(NeutGeneSummary$GeneName %in% Enriched$GeneName) & !(NeutGeneSummary$GeneName %in% Depleted$GeneName)) ,]

#Subsetting Data to Exclude sgRNAs that Were Discordant with the Category Assigned to Each Gene
EnrichedGenes<-Data_Ordered[which(Data_Ordered$GeneName %in% Enriched$GeneName),]
EnrichedGenes_UPsgRNAs<-EnrichedGenes[which(EnrichedGenes$sgRNA_Up==1),]
EnrichedGenes_UPsgRNAs$WeightedLogits<-(EnrichedGenes_UPsgRNAs$SortVsUnsort_NormRankDiff)*log(EnrichedGenes_UPsgRNAs$sgRNA_pvalue/(1-EnrichedGenes_UPsgRNAs$sgRNA_pvalue))
DepletedGenes<-Data_Ordered[which(Data_Ordered$GeneName %in% Depleted$GeneName),]
DepletedGenes_DOWNsgRNAs<-DepletedGenes[which(DepletedGenes$sgRNA_Down==1),]
DepletedGenes_DOWNsgRNAs$WeightedLogits<-(DepletedGenes_DOWNsgRNAs$SortVsUnsort_NormRankDiff)*log(DepletedGenes_DOWNsgRNAs$sgRNA_pvalue/(1-DepletedGenes_DOWNsgRNAs$sgRNA_pvalue))
NeutGenes<-Data_Ordered[which(Data_Ordered$GeneName %in% Neutral$GeneName),]
NeutGenes_NEUTsgRNAs<-NeutGenes[which(NeutGenes$sgRNA_Neutral==1),]
NeutGenes_NEUTsgRNAs$WeightedLogits<-(NeutGenes_NEUTsgRNAs$SortVsUnsort_NormRankDiff)*log(NeutGenes_NEUTsgRNAs$sgRNA_pvalue/(1-NeutGenes_NEUTsgRNAs$sgRNA_pvalue))
Data_Filtered<-rbind(EnrichedGenes_UPsgRNAs,NeutGenes_NEUTsgRNAs,DepletedGenes_DOWNsgRNAs)


#This step deploys the weighted logit method for combining p-values to arrive at a gene-level score by combining evidence across the sgRNAs targetting the 
#gene that are present after the filtering step above. Mudholkar and George (see https://apps.dtic.mil/dtic/tr/fulltext/u2/a049993.pdf) showed that the sum 
#of logit(p-value) takes the form of a scaled Student's t-distribution with 5m+4 degrees of freedom (df) and the scalar C=sqrt((0.3*df)/(m*(df-2))), where m  
#is the number of p-values combined. When including weights, the sum of Weight*logit(p) also follows a T-distribution, but df and C become a function of the 
#weights. Applying Mudholkar and George's closed form solution for C and df would therefore involve calculating a unique df and scaler at each gene. For 
#computational simplicity the implementation below instead assigns a conservative df=2 and uses an L2E estimate to determine an optimal C across the genes scored.
#We find df=2 to be conservative as the unweighted method would have df=14 for combining p-values from just two sgRNA and most screens aim for 3+ active sgRNAs
#of genes scored. When calculating our score, we weight logit(p) by the normalized change in rank to arrive at a combined score that accounts for significance  
#(p-value) and effect size (change in rank) at each sgRNA. The sum/C yields a gene level T-statisitic we've labeled SLIDER scores. p-values for each gene are 
#then computed from a Student's t-distribution with df=2. Multiple comparisons are corrected for using a Benjamini-Hochberg estimate of the FDR at each gene. 
#This script scores a single unsorted-vs-sorted replicate. To facilitate generating screen-level summary scores by combining replicates using Stouffer's method,
#a gene level Zequiv score computed by taking the inverse-normal(Gene_p-value) and matching its sign to the SLIDER-score.

GeneStat<-aggregate(Data_Filtered$WeightedLogits, by=list(Data_Filtered$GeneName),FUN=sum)
GeneStat$GeneName<-GeneStat$Group.1
GeneStat$SLIDER_NegWeightedLogitSum<-(-1*GeneStat$x)
GeneStat$Group.1<-NULL
GeneStat$x<-NULL
Gene_meanRankDiff<-aggregate(Data_Filtered$SortVsUnsort_NormRankDiff, by=list(Data_Filtered$GeneName),FUN=mean)
Gene_meanRankDiff$GeneName<-Gene_meanRankDiff$Group.1
Gene_meanRankDiff$meanDiffRank<-Gene_meanRankDiff$x
Gene_meanRankDiff$Group.1<-NULL
Gene_meanRankDiff$x<-NULL
GeneScores<-merge(x=GeneStat,y=Gene_meanRankDiff,by='GeneName')

#Using the 99% of Neutral Genes with Weighted Logit Sums Nearest to the Median to Inform an L2E Estimate of the 
#scalar of the Negative Weigthed Logit Sum that best approximates a t-distribution with 2-degrees of freedom.
#This is done to limit screen hits (not neutral) and outliers (i.e. neutral genes with read out issues) from 
#biasing the value estimated for the scalar.
L2E_neuts<-GeneScores$SLIDER_NegWeightedLogitSum[which(GeneScores$GeneName %in% NeutGenes$GeneName)]
L2E_neuts_ordered<-sort(L2E_neuts)
StartCount<-round((0.005*length(L2E_neuts_ordered)),0)
EndCount<-round((0.995*length(L2E_neuts_ordered)),0)
L2E_presubset<-L2E_neuts_ordered[StartCount:EndCount] - median(L2E_neuts_ordered[StartCount:EndCount])
if (abs(min(L2E_presubset))<abs(max(L2E_presubset))) {
  L2E_subset<-L2E_presubset[which(abs(L2E_presubset)<abs(min(L2E_presubset)))]
}
if (abs(min(L2E_presubset))>abs(max(L2E_presubset))) {
  L2E_subset<-L2E_presubset[which(abs(L2E_presubset)<abs(max(L2E_presubset)))]
}
degfree<-2
params<-c(0.2)
names(params) <- c('scalar')
obj=function(S){
  ((sqrt(pi)*gamma(x=((1/2)+(degfree))))/(((degfree)^(3/2))*S['scalar']*(beta(a=(degfree/2),b=(1/2)))^2*gamma(x = degfree)))-(2*mean((1/S['scalar'])*dt((L2E_subset)/S['scalar'],df=degfree)))
}
SLIDER_Tdist_params<-nlminb(params,obj, lower=0.01, upper=5)$par

#Calculating SLIDER_score T-statistic, determining p-value, correcting for multiple comparisons with FDR determined using the 
#Benjamini-Hochberg method, and reporting an inverse-normal(p-value) Zscore that would achieve the same p-value in a two-tail test. 
GeneScores$SLIDER_score<-(GeneScores$SLIDER_NegWeightedLogitSum-median(L2E_neuts_ordered[StartCount:EndCount]))/SLIDER_Tdist_params[1]
GeneScores$pvalue<-2*pt(q=(-1*abs(GeneScores$SLIDER_score)), df = degfree, lower.tail = TRUE, log.p=FALSE)
GeneScores$FDR<-p.adjust(GeneScores$pvalue,method = 'fdr')
GeneScores$Zequiv<-qnorm(GeneScores$pvalue/2)*(-1)*sign(GeneScores$SLIDER_score)

#Exporting Output Files. sgRNA level scoring and categroization are output into a seperate CSV file from gene level scores. 
filename<-paste(StudyName,'SLIDER_L2E_Tdist.txt')
sink(filename)
print(SLIDER_Tdist_params)
sink()
Filename<-paste(StudyName,'_SLIDER_sgRNA_scores_categorization.csv')
write.csv(x = Data_Ordered, Filename)
Filename<-paste(StudyName,'_SLIDER_GeneScores.csv')
write.csv(x = GeneScores, Filename)
Filename<-paste(StudyName,'_SLIDER_EnrichedGenes_ActiveGuides.csv')
write.csv(x = Enriched, Filename)
Filename<-paste(StudyName,'_SLIDER_DepletedGenes_ActiveGuides.csv')
write.csv(x = Depleted, Filename)
Filename<-paste(StudyName,'_SLIDER_NeutralGenes_NeutralGuides.csv')
write.csv(x = Neutral, Filename)

#Outputting Plots to a PDF file
Filename<-paste(StudyName,'_SLIDER_plots.pdf')
pdf(file = Filename)

#Plotting Distribution of Normalized Change in sgRNA Rank between the sorted and unsorted sample
plotname<-paste(StudyName,'\n Histogram of Normalized Change in sgRNA Rank')
b_graph <- sum(abs(Data$SortVsUnsort_NormRankDiff - 0))/length(Data$SortVsUnsort_NormRankDiff)
ymax_plot <- round((1/(2*b_graph)),1)
x_graph <- seq(-1,1,0.01)
laplace_line <-(1/(2*b_graph))*exp(-abs(x_graph)/b_graph)
hist(Data$SortVsUnsort_NormRankDiff, breaks=75, xlim=c(-1,1), col=rgb(0,1,0,0.4), ylim=c(0,ymax_plot), main=plotname, xlab='Normalized Change in sgRNA Rank', prob=T)
lines(x_graph, laplace_line, type='l', col='darkgreen', lwd=2)
legend('topright', c('Normalized Change in Rank', 'Laplace Distribution'),col=c('green','darkgreen'), bty='y',lty=c(1,1,1), cex=0.75)

#Plotting Laplace b parameter as a function of normalized rank in the unsorted sample (will appear as two lines, one for b+ and one for b-)
plotname<-paste(StudyName,'\n Laplace Parameter "b" vs Unsorted sgRNA Rank')
plot(Data_Ordered$Unsort_NormRank, Data_Ordered$b, xlab = "sgRNA Rank Prior to Sorting", ylab="Laplace Parameter b", main=plotname, pch=20)

#Plotting Distibution of sgRNA Abundance (Log-2-Normalized Read Counts) Before and After Sorting  
plotname<-paste(StudyName,'\n sgRNA Abundance Before and After Sorting')
hist(Data_Ordered$Unsort_log2PC, breaks=100, xlim=c(0,10), col=rgb(0,0,1,0.4), main=plotname, xlab='sgRNA Abundance (log2(PC))', prob=T, ylim=c(0,0.75))
hist(Data_Ordered$Sort_log2PC, breaks=100, xlim=c(0,10), col=rgb(1,0,0,0.4), prob=T, add=T)
legend('topright', c('Unsorted', 'Sorted'),col=c('blue','red'), bty='y',lty=c(1,1,1), cex=0.75)

#Plotting -log10(sgRNA_pvalue) Colormapped Scatter Plots and an Image Figure of Colormap
z=matrix(1:20,nrow=1)
x=1
y=c(0,0.75,1.5,max(-log10(Data_Ordered$sgRNA_pvalue)))
y=append(y,seq(from=2.0, to=5, len=17),after=3)
redColRamp <- colorRampPalette(c('black','red'))
sgRNApvalue_Col<- redColRamp(20)[as.numeric(cut(-log10(Data_Ordered$sgRNA_pvalue), breaks = y))]
plotname<-paste(StudyName,'\n Sorted vs Unsorted sgRNA Abundance')
plot(Data_Ordered$Unsort_log2PC, Data_Ordered$Sort_log2PC, col=sgRNApvalue_Col, xlab = "Unsorted (log2(PC))", ylab="Sorted (log2(PC))", main=plotname, pch=20)

#Plotting -log10(sgRNA_pvalue) Plots of Unsort-vs-Sort Rank and an Image Figure of Colormap
z=matrix(1:20,nrow=1)
x=1
y=c(0,0.75,1.5,max(-log10(Data_Ordered$sgRNA_pvalue)))
y=append(y,seq(from=2.0, to=5, len=17),after=3)
redColRamp <- colorRampPalette(c('black','red'))
sgRNApvalue_Col<- redColRamp(20)[as.numeric(cut(-log10(Data_Ordered$sgRNA_pvalue), breaks = y))]
plotname<-paste(StudyName,'\n Sorted vs Unsorted sgRNA Ranks')
plot(Data_Ordered$Unsort_NormRank, Data_Ordered$Sort_NormRank, col=sgRNApvalue_Col, xlab = "Unsorted sgRNA Rank", ylab="Sorted sgRNA Rank", main=plotname, pch=20)

#plotting LFC against Unsorted Rank for each sgRNA as a -log10(sgRNA_pvalue) colormapped scatter plot
z=matrix(1:20,nrow=1)
x=1
y=c(0,0.75,1.5,max(-log10(Data_Ordered$sgRNA_pvalue)))
y=append(y,seq(from=2.0, to=5, len=17),after=3)
redColRamp <- colorRampPalette(c('black','red'))
sgRNApvalue_Col<- redColRamp(20)[as.numeric(cut(-log10(Data_Ordered$sgRNA_pvalue), breaks = y))]
plotname<-paste(StudyName,'\n Unsorted sgRNA Rank vs Log Fold Change')
plot(Data_Ordered$Unsort_NormRank, (Data_Ordered$Sort_log2PC - Data_Ordered$Unsort_log2PC), col=sgRNApvalue_Col, xlab = "Unsorted sgRNA Rank", ylab="sgRNA LFC Sorted-vs-Unsorted", main=plotname, pch=20)

#plotting LFC against Sorted Rank for each sgRNA as a -log10(sgRNA_pvalue) colormapped scatter plot
z=matrix(1:20,nrow=1)
x=1
y=c(0,0.75,1.5,max(-log10(Data_Ordered$sgRNA_pvalue)))
y=append(y,seq(from=2.0, to=5, len=17),after=3)
redColRamp <- colorRampPalette(c('black','red'))
sgRNApvalue_Col<- redColRamp(20)[as.numeric(cut(-log10(Data_Ordered$sgRNA_pvalue), breaks = y))]
plotname<-paste(StudyName,'\n Sorted sgRNA Rank vs Log Fold Change')
plot(Data_Ordered$Sort_NormRank, (Data_Ordered$Sort_log2PC - Data_Ordered$Unsort_log2PC), col=sgRNApvalue_Col, xlab = "Sorted sgRNA Rank", ylab="sgRNA LFC Sorted-vs-Unsorted", main=plotname, pch=20)

#plotting image of colormap for sgRNA pvalues
y2<-c(0,0.75,1.5,seq(from=2.0, to=5, len=17),6)
image(x,y2, z,col=redColRamp(20), ylab="-log10(sgRNA_pval)")

#Plotting Each Gene's SLIDER Score Against the Mean Difference in Rank of its sgRNAs When Comparing Sorted to Unsorted
z=matrix(1:20,nrow=1)
x=1
y=c(0,0.75,1.5,(max(-log10(GeneScores$pvalue))+1))
y=append(y,seq(from = 2.0, to = (max(-log10(GeneScores$pvalue))-1),len=17),after=3)
redColRamp <- colorRampPalette(c('black','red'))
Gene_pvalue_Col<- redColRamp(20)[as.numeric(cut(-log10(GeneScores$pvalue), breaks = y))]
plotname<-paste(StudyName,'\n SLIDER_Tstat vs Average Change in Rank of sgRNA for Each Gene Colored by pvalue')
plot(GeneScores$meanDiffRank, GeneScores$SLIDER_score, xlim=c(-0.6,0.6), ylim=c(-70,70), col=Gene_pvalue_Col, xlab = "Mean Difference in Rank (Sorted-vs-Unsorted)", ylab="SLIDER Score", main=plotname, pch=20)

#plotting colormap of Gene pvalue
image(x,y, z,col=redColRamp(20), ylab = "-log(Gene_pval)")

#Plotting Each Gene's -log10(pvalue) Against the Mean Change in Rank of its sgRNA from Comparing Sorted to Unsorted (Rank Based Volcano Plot)
z=matrix(1:20,nrow=1)
x=1
y1=seq(from = -35, to = 35,len=18)
minval<-(min(GeneScores$SLIDER_score)-1)
y=append(minval,y1,after=length(minval))
y=append(y,(max(GeneScores$SLIDER_score)+1),after=length(y))
redColRamp <- colorRampPalette(c('blue','black','red'))
SLIDER_T_Col<- redColRamp(20)[as.numeric(cut(GeneScores$SLIDER_score, breaks = y))]
plotname<-paste(StudyName,'\n Volcano Plot for Genes Assayed (by Rank), Colored by SLIDER Score')
plot(GeneScores$meanDiffRank, -log10(GeneScores$pvalue), col=SLIDER_T_Col, xlab = "Mean Change in Rank (Sorted-vs-Unsorted)", ylab="-log10(pvalue)", main=plotname, pch=20, xlim=c(-0.75,0.75), ylim=c(0,(1+(-log10(min(GeneScores$pvalue))))))

#Plotting colormap of SLIDER stats
y2<-append((min(y1)-5),y1,after=length(min(y1)-3))
y2<-append(y2,max(y2+5),after=length(y))
image(x,y2,z,col=redColRamp(20), ylab = "SLIDER score")

#Plotting Tdist Over Histogram of SLIDER Scores 
plotname<-paste(StudyName,'\n t-Distirbution Over Histogram of SLIDER Scores')
hist(GeneScores$SLIDER_score, breaks=2000, xlim=c(-10,10), col=rgb(0,1,0,0.4), main=plotname, xlab='Statistic', prob=T, ylim=c(0,0.7))
curve(dt(x, df=degfree), col='darkgreen', lwd=2, add=TRUE, yaxt='n', n=1001)
legend('topright', c('SLIDER Score', 'Students t-Dist'),col=c('green','darkgreen'), bty='y',lty=c(1,1,1), cex=0.75)

#Generating Dot Plots of SLIDER Scores Ranked from Largest to Smallest and an Image of the Matching Colormap. Making Two plots, one with and one without gene name labels.
plotname<-paste(StudyName,'\n Ranked SLIDER Scores')
RankedGeneScores<-GeneScores[order(-GeneScores$SLIDER_score),]
z=matrix(1:20,nrow=1)
x=1
y=c(0,0.75,1.5,(max(-log10(RankedGeneScores$pvalue))+1))
y=append(y,seq(from = 2.0, to = (max(-log10(RankedGeneScores$pvalue))-1),len=17),after=3)
redColRamp <- colorRampPalette(c('black','red'))
RankedGene_pvalue_Col<- redColRamp(20)[as.numeric(cut(-log10(RankedGeneScores$pvalue), breaks = y))]
plot(RankedGeneScores$SLIDER_score, col=RankedGene_pvalue_Col, xlab="Rank by SLIDER score", ylab="SLIDER Score", main=plotname, pch=20, xlim=c(0,(2000+length(RankedGeneScores$SLIDER_score))))
plot(RankedGeneScores$SLIDER_score, col=RankedGene_pvalue_Col, xlab="Rank by SLIDER score", ylab="SLIDER Score", main=plotname, pch=20, xlim=c(0,(2000+length(RankedGeneScores$SLIDER_score))))
text(Indices<-seq(from=1, to=length(RankedGeneScores$SLIDER_score), by = 1), RankedGeneScores$SLIDER_score, labels = RankedGeneScores$GeneName, cex=0.6, pos=4)
image(x,y,z,col=redColRamp(20), ylab = "Gene -log10(pvalue)")

dev.off()
