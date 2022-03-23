rm(list=ls())
### Pipeline adapted from DADA2 tutorial by Varun Gupta
### make sure you have txt file with list of sample ### names used for files

##### DADA pipeline  #####

## downloading dada2 package

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")
BiocManager::install(c("dada2"))
library(dada2)
library(Rcpp)
file.choose()
path <- file.choose
list.files(path)

## Now we read in the names of the fastq files, and perform some string manipulation to get matched lists of the forward and reverse fastq files.##
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#Inspect read quality profiles
##We start by visualizing the quality profiles of the forward and reverse reads:
plotQualityProfile(fnFs[1:7])
#The forward read quality reduces after 240 bp. So will need to cut after that ##
plotQualityProfile(fnRs[1:7])
#The reserve read quality reduces after 160 bp. So will need to cut after that ##

####Filter and trim
#Assign the filenames for the filtered fastq.gz files.#
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
#  maxN should be set to 0. It removes any ambigous bp. dada pipeline cannot handle them#
#  maxEE: maximum number of errors in a read. 2 is a good baseline. If the reads are of poor quality, then should consider makinh the number higher#
#  truncQ: truncates reads at the first instance of a Q score less than or equal to the value specified. Q of 2 is very bad score. Do not lower it#
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
head(out)

####  Learn the Error Rates  ####
#Every batch will have a different error rate #
# First run on the forward reads #
errF <- learnErrors(filtFs, multithread=TRUE)
# Then run on the reverse reads #
errR <- learnErrors(filtRs, multithread=TRUE)
#It is always worthwhile to visualize the estimated error rates:#
# See if there is a goot fit between the black (actual error rate) and red line (error rate calculated with Q score)#
plotErrors(errF, nominalQ=TRUE)

#####   Dereplication    ###
#Dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance” equal to the number of reads with that unique sequence. Dereplication substantially reduces computation time by eliminating redundant comparisons.#
#DADA2 retains a summary of the quality information associated with each unique sequence.#
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#######    Sample Inference    ######
# Using the error model developed earlier, the algorithm calculates abundance p values for each unique sequence#
#Tests null hypothesis that a sequence with a given error rate is too abundant to be explained by sequencing error#
#Low p-value: actually a read read, not caused by random errors#
#High p-value: a read was likely caused by errors#
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
#Inspecting the returned dada-class object#
dadaFs[[1]]
dadaRs[[1]]

#######     Merge Paired End Reads    #####
## Will merge if forward and reverse reads overlap by atleast 12 bases ##
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

####### Construct sequence table    #####
#Construct a ASV table #
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

###### Remove Chimeras   #####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
# This will tell what percentage of your total reads are not a chimera #
sum(seqtab.nochim)/sum(seqtab)

##### Track reads through the pipeline  #####
# Checks for the number of reads that made it through each step in the pipeline #
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

##### Assign Taxonomy   ######
# I downloaded Silva 132 database#
# Also downloaded Silva 132 species level assignment #
# This command only goes down till Genus level #
taxa <- assignTaxonomy(seqtab.nochim, "~/Desktop/Megan/silva_nr_v138_train_set.fa.gz", multithread=TRUE)
#Let’s inspect the taxonomic assignments #
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#####phylogenetic tree#####
> retrieveseqs <- function(seqnames,acnucdb)
{
  myseqs <- list()   # Make a list to store the sequences
  require("seqinr")  # This function requires the SeqinR R package
  choosebank(acnucdb)
  for (i in 1:length(seqnames))
  {
    seqname <- seqnames[i]
    print(paste("Retrieving sequence",seqname,"..."))
    queryname <- "query2"
    query <- paste("AC=",seqname,sep="")
    query(`queryname`,`query`)
    seq <- getSequence(query2$req[[1]]) # Makes a vector "seq" containing the sequence
    myseqs[[i]] <- seq
  }
  closebank()
  return(myseqs)
}



## R package 3.5.0 or greater needed for dada2 pipeline ##
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")


library(dada2)

path <- "~/Desktop/AS.DATA"
list.files(path)
### Now we read in the names of the fastq files, and perform some string manipulation to get matched lists of the forward and reverse fastq files.##
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

####    Inspect read quality profiles     #####
##We start by visualizing the quality profiles of the forward and reverse reads:
plotQualityProfile(fnFs[1:2])
#The forward read quality reduces after 200 bp. So will need to cut after that ## <10,000 reads is very low for Illumina, consider re-sequencing
plotQualityProfile(fnRs[1:2])
#The reserve read quality reduces after 190 bp. So will need to cut after that ##

####     Filter and trim     ####
#Assign the filenames for the filtered fastq.gz files.#
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
#  maxN should be set to 0. It removes any ambigous bp. dada pipeline cannot handle them#
#  maxEE: maximum number of errors in a read. 2 is a good baseline. If the reads are of poor quality, then should consider making the number higher#
#  truncQ: truncates reads at the first instance of a Q score less than or equal to the value specified. Q of 2 is very bad score. Do not lower it#
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(180,180),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
head(out)

####  Learn the Error Rates  ####
#Every batch will have a different error rate #
# First run on the forward reads #
errF <- learnErrors(filtFs, multithread=TRUE)
# Then run on the reverse reads #
errR <- learnErrors(filtRs, multithread=TRUE)
##It is always worthwhile to visualize the estimated error rates:#
##Note: In a large data set, the loop will stop after 10 cycles. Check if see if the the below conditions are met, if not, you can force the function to do more cycles#
# Check if the last two values are much much smaller than the first values # and make sure that the values at the end approach zero
dada2:::checkConvergence(errF)
dada2:::checkConvergence(errR)
# See if there is a goot fit between the black dot (actual error rate) and black line ( fitted error rate)#
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#####   Dereplication    ###
#Dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance” equal to the number of reads with that unique sequence. Dereplication substantially reduces computation time by eliminating redundant comparisons.#
#DADA2 retains a summary of the quality information associated with each unique sequence.#
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#######    Sample Inference    ######
# Using the error model developed earlier, the algorithm calculates abundance p values for each unique sequence#
#Tests null hypothesis that a sequence with a given error rate is too abundant to be explained by sequencing error#
#Low p-value: actually a read read, not caused by random errors#
#High p-value: a read was likely caused by errors#
#Should be clustering down to fewer unique sequences, in better data#
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool = "pseudo")
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool = "pseudo")
#Inspecting the returned dada-class object# seeing how many sequence variances there are & the forward and reverse should be very similar
dadaFs[[1]]
dadaRs[[1]]

#######     Merge Paired End Reads    #####
## Will merge if forward and reverse reads overlap by atleast 12 bases ##\
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

####### Construct sequence table    #####
#Construct a ASV table # 
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

###### Remove Chimeras   #####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
# This will tell what percentage of your total reads are not a chimera - if losing over 10-15%, data is bad #
# Chimera = an organism or tissue that contains at least two different sets of DNA, most often originating from the fusion of as many different zygotes (fertilized eggs) #
sum(seqtab.nochim)/sum(seqtab)

##### Track reads through the pipeline  #####
# Checks for the number of reads that made it through each step in the pipeline #
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

##### Assign Taxonomy   ######

# I downloaded Silva 132 database#
# Also downloaded Silva 132 species level assignment #
# This command only goes down till Genus level # 
# Bootstrap value has to equal 50, in order to assign taxonomy, can however raise or lower that bootstrap value #
taxa <- assignTaxonomy(seqtab.nochim, file.choose(), multithread=TRUE)
# This adds species to taxanomic assignment to above #
taxaplus <- addSpecies(taxa, "~/Desktop/silva_species_assignment_v132.fa", verbose=TRUE)
#the silva databases have been updated since this 
####### End of DADA2 pipeline  #######



library(phyloseq)
##We now construct a phyloseq object directly from the dada2 outputs.
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa))
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps

sample_data(ps)$samples.out = factor(as.character(sample_data(ps)$samples.out), levels = c("40E", "40P", "50E", "50S", "BLNK"))
test<-otu_table(ps,taxa_are_rows(samples.out))


test1<-data.frame(test)

write.csv(tax_table(ps), "taxtable.csv")
write.csv(otu_table(ps), "otutable.csv")
#### Handoff to Phyloseq and ggplot2 #####

library(phyloseq)
library(ggplot2)
library(digest)
library(dplyr)
library(RColorBrewer)
library(readr)

## Importing mapping file into phyloseq  ####
library(readr)
Antarctic_Endoliths_sample_metadata_2 <- read_csv("~/Desktop/LL/Antarctic Endoliths_sample metadata-2.csv", 
                                                  col_types = cols(X5 = col_skip(), X6 = col_skip()))
View(Antarctic_Endoliths_sample_2)
map <- (Antarctic_Endoliths_sample_2)

#map$Year <- as.factor(as.character(map$Year))
### get sample IDs to match (strip off sequencer sample # from mapping file) ####

map$S <- as.character(map$Site)
samplesite <- strsplit(map$Site, split = "\\.S")
for(i in 1:length(samplesite)){
  map$Site[i] <- samplesite[[i]][1]
}

sample_names(map) <- map$Name

map <- sample_data(map)
samplenames(map) <- map$Name
sample_names(map)

## get relative abundance
ps <- transform_sample_counts(ps, function(x) x/sum(x))
ps = filter_taxa(ps, function(x) var(x) > 1e-05, TRUE)

## same as above but >1% OTUs ##
ps<- filter_taxa(ps, function(OTU) sum(OTU) > 0.01, TRUE)
ps <- transform_sample_counts(ps, function(OTU) 100* OTU/sum(OTU)) 


####### End of DADA2 pipeline  #######
#### Handoff to Phyloseq and ggplot2 #####
library(dplyr)
library(magrittr)
library(pillar)
library(genomes)
library(readr)
library(curl)
library(GenomicRanges)
library(GenomeInfoDb)
library(ggpubr)
library(dada2)
library(Rcpp)
library(phyloseq)







####<1% abundance (Make 100% bar plots with less abundance in a sub-category)#### 
glomk2 <- tax_glom(ps,taxrank = 'Kingdom')
data_glomk2 <- psmelt(glomk2)
data_glomk2$Kingdom <- as.character(data_glomk2$Kingdom)
data_glomk2$Kingdom[data_glomk2$Abundance < 0.01] <- "< 0.5% Abundance"

glomp2 <- tax_glom(ps,taxrank = 'Phylum')
data_glomp2 <- psmelt(glomp2)
data_glomp2$Phylum <- as.character(data_glomp2$Phylum)
data_glomp2$Phylum[data_glomp2$Abundance < 0.01] <- "< 1% Abundance"

glomc2 <- tax_glom(ps.site,taxrank = 'Class')
data_glomc2 <- psmelt(glomc2)
data_glomc2$Class <- as.character(data_glomc2$Class)
data_glomc2$Class[data_glomc2$Abundance < 0.01] <- "< 1% Abundance"

glomo2 <- tax_glom(ps.w.site,taxrank = 'Order')
data_glomo2 <- psmelt(glomo2)
data_glomo2$Order <- as.character(data_glomo2$Order)
data_glomo2$Order[data_glomo2$Abundance < 0.01] <- "< 1% Abundance"

glomf2 <- tax_glom(ps.w.site,taxrank = 'Family')
data_glomf2 <- psmelt(glomf2)
data_glomf2$Family <- as.character(data_glomf2$Family)
data_glomf2$Family[data_glomf2$Abundance < 0.01] <- "< 1% Abundance"

glomg2 <- tax_glom(ps.w.site,taxrank = 'Genus')
data_glomg2 <- psmelt(glomg2)
data_glomg2$Genus <- as.character(data_glomg2$Genus)
data_glomg2$Genus[data_glomg2$Abundance < 0.01] <- "< 1% Abundance"

gloms2 <- tax_glom(ps.w.site,taxrank = 'Species')
data_gloms2 <- psmelt(gloms2)
data_gloms2$Species <- as.character(data_gloms2$Species)
data_gloms2$Species[data_gloms2$Abundance < 0.01] <- "< 1% Abundance"




####MERGE WITH MAPPING FILE####
## Do this before merging with a phyloseq file ##
map$S <- as.character(map$S)
Ss <- strsplit(map$S, split = "\\.S")
for(i in 1:length(Ss)){
  map$S[i] <- Ss[[i]][1]
}

data_glomk <- merge.data.frame(data_glomk, map)
data_glomk

data_glomp <- merge.data.frame(data_glomp2, map)
data_glomp

data_glomc <- merge.data.frame(data_glomc2, map)
data_glomc

data_glomo <- merge.data.frame(data_glomo, map)
data_glomo

data_glomf <- merge.data.frame(data_glomf, map)
data_glomf

data_glomg <- merge.data.frame(data_glomg, map)
data_glomg

######FINAL CODE PLOTS######

MyKingdomPalette <- c("darkorchid2", "orangered", "springgreen", "pink1", "darkorchid2")
KingdomGraph <- ggplot(data=data_glomk, aes(x=Sample, y=Abundance, fill=Kingdom))+
  geom_bar(stat="identity", position="stack")+ scale_fill_manual(values = MyKingdomPalette) +
  theme_bw() + labs(x = "Sample Name", y = "Relative Abundance (%)") +theme(axis.text.x = element_text(angle = 90))
print(KingdomGraph)


MyPhylumPalette <- c("grey", "blue1", "lightblue", "orangered", "springgreen", "pink1", "darkorchid2", "turquoise", "red", "pink3", "violetred1", "orange", "lightpink", "blue", "peachpuff", "black", "darkred", "white", "plum", "green")
PhylumGraph <- ggplot(data=data_glomp, aes(x=Sample, y= Abundance, fill=Phylum))+ 
  geom_bar(stat="identity", position="stack")+ scale_fill_manual(values = MyPhylumPalette) +
  theme_bw() + labs(x = "Sample Name", y = "Relative Abundance (%)") +theme(axis.text.x = element_text(angle = 90))
print(PhylumGraph)

MyClassPalette <- c("grey", "tomato", "darkmagenta", "yellow","lightblue", "orangered", "springgreen", "pink1", "darkorchid2", "turquoise", "red", "pink3", "violetred1", "orange", "lightpink", "blue", "peachpuff", "purple", "springgreen1", "aquamarine1", "darkorchid2", "orange", "mediumblue", "goldenrod1", "violetred1", "cornflowerblue", "dodgerblue", "red", "blue", "chocolate", "green", "violetred1", "black", "gold", "violet", "darkslateblue", "cadetblue1", "olivedrab", "palegoldenrod", "grey")
ClassGraph <- ggplot(data=data_glomc2, aes(x=Sample, y=Abundance, fill=Class))+
  geom_bar(stat="identity", position="stack")+ scale_fill_manual(values = MyClassPalette) +
  theme_bw() + labs(x = "Sample Name", y = "Relative Abundance (%)") +theme(axis.text.x = element_text(angle = 90))
print(ClassGraph)


#Order 1% combination#
MyOrderPalette <- c("grey", "tomato", "darkmagenta", "yellow", "lightblue", "orangered", "springgreen", "pink1", "darkorchid2", "turquoise", "red", "pink3", "violetred1", "darkgreen", "lightpink", "blue", "peachpuff", "purple", "springgreen1", "orchid", "darkorchid2", "firebrick4", "mediumblue", "goldenrod1", "violetred1", "cornflowerblue", "dodgerblue", "red", "blue", "chocolate", "green", "violetred1", "yellow", "gold", "violet","slateblue1", "tomato", "darkmagenta", "purple", "goldenrod1", "lemonchiffon1", "springgreen1", "maroon", "maroon2", "mediumblue", "gold", "dodgerblue", "orange", "mistyrose", "darkorange", "purple4", "red2", "greenyellow", "blue", "peachpuff", "goldenrod1", "violetred1", "cornflowerblue", "dodgerblue", "red", "blue", "chocolate", "peachpuff", "violetred1", "greenyellow", "gold", "violet", "darkslateblue", "black", "olivedrab", "palegoldenrod", "darkgreen", "steelblue1", "goldenrod1", "violetred1", "cornflowerblue", "dodgerblue", "red", "blue", "chocolate", "peachpuff", "violetred1", "greenyellow", "gold", "violet", "darkslateblue", "black", "olivedrab", "palegoldenrod","purple", "yellow3", "plum", "palegreen2", "lightcoral", "chartreuse1", "maroon", "slateblue1", "tomato", "purple", "darkred", "hotpink1", "skyblue1", "darkgreen", "lemonchiffon1", "springgreen1", "maroon2", "dodgerblue", "mediumblue", "yellow3", "coral", "darkolivegreen4", "midnightblue", "orangered", "slateblue", "khaki","seagreen" , "mediumpurple2", "white", "greenyellow", "darkred", "olivedrab", "yellow", "forestgreen", "grey", "violet")
OrderGraph <- ggplot(data=data_glomo, aes(x=Sample, y=Abundance, fill=Order))+
  geom_bar(stat="identity", position="stack")+ scale_fill_manual(values = MyOrderPalette) +
  theme_bw() + labs(x = "Sample Name", y = "Relative Abundance (%)") +theme(axis.text.x = element_text(angle = 90))
print(OrderGraph)

#Family 1% combo#
MyFamilyPalette <- c("grey", "darkgreen", "steelblue1", "goldenrod1", "violetred1", "cornflowerblue", "dodgerblue", "red", "blue", "chocolate", "peachpuff", "violetred1", "greenyellow", "gold", "violet", "darkslateblue", "cadetblue1", "olivedrab", "palegoldenrod","purple", "yellow3", "plum", "palegreen2", "lightcoral", "chartreuse1", "maroon", "slateblue1", "tomato", "purple", "darkred", "hotpink1", "skyblue1", "darkgreen", "lemonchiffon1", "springgreen1", "maroon2", "dodgerblue", "mediumblue", "yellow3", "coral", "darkolivegreen4", "chartreuse1", "slateblue1", "tomato", "darkmagenta", "purple", "goldenrod1", "lemonchiffon1", "chartreuse1", "maroon", "slateblue1", "tomato", "purple", "darkred", "hotpink1", "skyblue1", "darkgreen", "lemonchiffon1", "springgreen1", "maroon2", "DODGERBLUE", "mediumblue", "yellow3", "coral", "darkolivegreen4", "cornflowerblue", "darkslategrey", "seagreen3", "orange", "darkorchid1", "mistyrose", "cadetblue1", "purple4", "red", "green4", "blue", "plum", "peachpuff", "midnightblue", "indianred3", "khaki", "mediumpurple2", "steelblue1", "greenyellow", "deepskyblue2", "olivedrab", "yellow", "black", "forestgreen", "darkslateblue", "rosybrown1", "grey", "navy", "palegreen")
FamilyGraph <- ggplot(data=data_glomf, aes(x=Sample, y=Abundance, fill=Family))+
  geom_bar(stat="identity", position="stack")+ scale_fill_manual(values = MyFamilyPalette) +
  theme_bw() + labs(x = "Sample Name", y = "Relative Abundance (%)") +theme(axis.text.x = element_text(angle = 90))
print(FamilyGraph)


#Genus <2% combo# - mistyrose1 (removed after steelblue1), sienna1 removed before navy
MyGenusPalette <- c("grey", "darkgreen", "steelblue1", "purple", "yellow3", "plum", "palegreen2", "lightcoral", "chartreuse1", "goldenrod1", "violetred1", "cornflowerblue", "dodgerblue", "red", "blue", "chocolate", "peachpuff", "violetred1", "greenyellow", "gold","maroon", "slateblue1", "tomato", "purple", "darkred", "hotpink1", "skyblue1", "darkgreen", "lemonchiffon1", "springgreen1", "maroon2", "dodgerblue", "mediumblue", "yellow3", "coral", "darkolivegreen4", "cornflowerblue", "darkslategrey", "seagreen3", "orange", "darkorchid1", "mistyrose", "cadetblue1", "purple4", "red", "green4", "blue", "plum", "peachpuff", "midnightblue", "indianred3", "khaki", "mediumpurple2", "steelblue1", "greenyellow", "deepskyblue2", "firebrick", "yellow", "lightcoral", "forestgreen", "darkslateblue", "rosybrown1", "grey", "navy", "palegreen", "mistyrose1", "slateblue1", "darkorchid1", "gold", "cadetblue1", "springgreen1", "maroon2", "dodgerblue", "rosybrown1", "darkolivegreen4", "skyblue1", "lavender", "orange","hotpink1", "GREENYELLOW", "LEMONCHIFFON1", "cyan1", "purple4", "red", "green4", "blue", "cornflowerblue", "olivedrab", "mediumpurple1", "darkred", "sandybrown", "navy", "violet", "khaki", "lawngreen", "paleturquoise1", "tan", "forestgreen", "yellow", "black", "tomato", "black", "grey", "sienna1")
GenusGraph <- ggplot(data=data_glomg, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(stat="identity", position="stack")+ 
  scale_fill_manual(values = MyGenusPalette) + 
  theme_bw() + labs(x = "Sample Name", y = "Relative Abundance (%)") +theme(axis.text.x = element_text(angle = 90))
print(GenusGraph)



#####COLS GROUPED BY SITE#####

MyKingdomPalette <- c("darkorchid2", "orangered", "springgreen", "pink1", "darkorchid2")
KingdomGraph <- ggplot(data= data_glomk2, aes(x=Sample, y=Abundance, fill=Kingdom))+
  geom_bar(stat="identity")+ scale_fill_manual(values = MyKingdomPalette) +
  theme_bw() + labs(x = "Sample Name", y = "Relative Abundance (%)") +theme(axis.text.x = element_text(angle = 90)) +
  scale_x_discrete(limits=c("3-S58", "37-S92", "2-S57", "4-S59", "25-S80", 
                            "5-S60", "6-S61", "17-S72", "33-S88", 
                            "9-S64", "14-S69", "28-S83", "24-S79", "16-S71", "19-S74", "30-S85", 
                            "1-S56", "15-S70", "20-S75", "23-S78", "22-S77", "21-S76", "26-S81", "31-S86", "27-S82",  "13-S68", "36-S91", 
                            "10-S65", "7-S62", "8-S63", "11-S66", "32-S87", 
                            "12-S67", "18-S73", "29-S84", "34-S89", "35-S90"))
print(KingdomGraph)




data_glomk <- data_glomk(sites = c("3-S58", "37-S92", "2-S57", "4-S59", "25-S80", 
                              "5-S60", "6-S61", "17-S72", "33-S88", 
                              "9-S64", "14-S69", "28-S83", "24-S79", "16-S71", "19-S74", "30-S85", 
                              "1-S56", "15-S70", "20-S75", "23-S78", "22-S77", "21-S76", "26-S81", "31-S86", "27-S82",  "13-S68", "36-S91", 
                              "10-S65", "7-S62", "8-S63", "11-S66", "32-S87", 
                              "12-S67", "18-S73", "29-S84", "34-S89", "35-S90")

library(ggplot2)
toplot.N <- data.frame(set=c("3-S58", "37-S92", "2-S57", "4-S59", "25-S80", 
                              "5-S60", "6-S61", "17-S72", "33-S88", 
                              "9-S64", "14-S69", "28-S83", "24-S79", "16-S71", "19-S74", "30-S85", 
                              "1-S56", "15-S70", "20-S75", "23-S78", "22-S77", "21-S76", "26-S81", "31-S86", "27-S82",  "13-S68", "36-S91", 
                              "10-S65", "7-S62", "8-S63", "11-S66", "32-S87", 
                              "12-S67", "18-S73", "29-S84", "34-S89", "35-S90"))
3-S58, 37-S92, 2-S57, 4-S59, 25-S80, 
5-S60, 6-S61, 17-S72, 33-S88, 
9-S64, 14-S69, 28-S83, 24-S79, 16-S71, 19-S74, 30-S85, 
1-S56, 15-S70, 20-S75, 23-S78, 22-S77, 21-S76, 26-S81, 31-S86, 27-S82, 13-S68, 36-S91, 
10-S65, 7-S62, 8-S63, 11-S66, 32-S87, 
12-S67, 18-S73, 29-S84, 34-S89, 35-S90))  
       


scale_x_discrete(limits=c(3-S58, 37-S92, 2-S57, 4-S59, 25-S80, 
                          5-S60, 6-S61, 17-S72, 33-S88, 
                          9-S64, 14-S69, 28-S83, 24-S79, 16-S71, 19-S74, 30-S85, 
                          1-S56, 15-S70, 20-S75, 23-S78, 22-S77, 21-S76, 26-S81, 31-S86, 27-S82, 13-S68, 36-S91, 
                          10-S65, 7-S62, 8-S63, 11-S66, 32-S87, 
                          12-S67, 18-S73, 29-S84, 34-S89, 35-S90))
########

MyPhylumPalette <- c("grey", "blue1", "lightblue", "orangered", "springgreen", "pink1", "darkorchid2", "turquoise", "red", "pink3", "violetred1", "orange", "lightpink", "blue", "peachpuff", "black", "darkred", "white", "plum", "green")
PhylumGraph <- ggplot(data=data_glomp, aes(x=Sample, y= Abundance, fill=Phylum))+ 
  geom_bar(stat="identity", position="stack")+ scale_fill_manual(values = MyPhylumPalette) +
  theme_bw() + labs(x = "Sample Name", y = "Relative Abundance (%)") +theme(axis.text.x = element_text(angle = 90))+
  scale_x_discrete(limits=c("3-S58", "37-S92", "2-S57", "4-S59", "25-S80", 
                            "5-S60", "6-S61", "17-S72", "33-S88", 
                            "9-S64", "14-S69", "28-S83", "24-S79", "16-S71", "19-S74", "30-S85", 
                            "1-S56", "15-S70", "20-S75", "23-S78", "22-S77", "21-S76", "26-S81", "31-S86", "27-S82",  "13-S68", "36-S91", 
                            "10-S65", "7-S62", "8-S63", "11-S66", "32-S87", 
                            "12-S67", "18-S73", "29-S84", "34-S89", "35-S90"))
print(PhylumGraph)

MyClassPalette <- c("grey", "tomato", "darkmagenta", "yellow","lightblue", "orangered", "springgreen", "pink1", "darkorchid2", "turquoise", "red", "pink3", "violetred1", "orange", "lightpink", "blue", "peachpuff", "purple", "springgreen1", "aquamarine1", "darkorchid2", "orange", "mediumblue", "goldenrod1", "violetred1", "cornflowerblue", "dodgerblue", "red", "blue", "chocolate", "green", "violetred1", "black", "gold", "violet", "darkslateblue", "cadetblue1", "olivedrab", "palegoldenrod", "grey")
ClassGraph <- ggplot(data=data_glomc2, aes(x=Sample, y=Abundance, fill=Class))+
  geom_bar(stat="identity", position="stack")+ scale_fill_manual(values = MyClassPalette) +
  theme_bw() + labs(x = "Sample Name", y = "Relative Abundance (%)") +
  theme(axis.text.x = element_text(angle = 90))+scale_x_discrete(limits=c("3-S58", "37-S92", "2-S57", "4-S59", "25-S80", 
                                                                                                                                              "5-S60", "6-S61", "17-S72", "33-S88", 
                                                                                                                                              "9-S64", "14-S69", "28-S83", "24-S79", "16-S71", "19-S74", "30-S85", 
                                                                                                                                              "1-S56", "15-S70", "20-S75", "23-S78", "22-S77", "21-S76", "26-S81", "31-S86", "27-S82",  "13-S68", "36-S91", 
                                                                                                                                              "10-S65", "7-S62", "8-S63", "11-S66", "32-S87", 
                                                                                                                                              "12-S67", "18-S73", "29-S84", "34-S89", "35-S90"))
print(ClassGraph)


#Order 1% combination#
MyOrderPalette <- c("grey", "tomato", "darkmagenta", "yellow", "lightblue", "orangered", "springgreen", "pink1", "darkorchid2", "turquoise", "red", "pink3", "violetred1", "darkgreen", "lightpink", "blue", "peachpuff", "purple", "springgreen1", "orchid", "darkorchid2", "firebrick4", "mediumblue", "goldenrod1", "violetred1", "cornflowerblue", "dodgerblue", "red", "blue", "chocolate", "green", "violetred1", "yellow", "gold", "violet","slateblue1", "tomato", "darkmagenta", "purple", "goldenrod1", "lemonchiffon1", "springgreen1", "maroon", "maroon2", "mediumblue", "gold", "dodgerblue", "orange", "mistyrose", "darkorange", "purple4", "red2", "greenyellow", "blue", "peachpuff", "goldenrod1", "violetred1", "cornflowerblue", "dodgerblue", "red", "blue", "chocolate", "peachpuff", "violetred1", "greenyellow", "gold", "violet", "darkslateblue", "black", "olivedrab", "palegoldenrod", "darkgreen", "steelblue1", "goldenrod1", "violetred1", "cornflowerblue", "dodgerblue", "red", "blue", "chocolate", "peachpuff", "violetred1", "greenyellow", "gold", "violet", "darkslateblue", "black", "olivedrab", "palegoldenrod","purple", "yellow3", "plum", "palegreen2", "lightcoral", "chartreuse1", "maroon", "slateblue1", "tomato", "purple", "darkred", "hotpink1", "skyblue1", "darkgreen", "lemonchiffon1", "springgreen1", "maroon2", "dodgerblue", "mediumblue", "yellow3", "coral", "darkolivegreen4", "midnightblue", "orangered", "slateblue", "khaki","seagreen" , "mediumpurple2", "white", "greenyellow", "darkred", "olivedrab", "yellow", "forestgreen", "grey", "violet")
OrderGraph <- ggplot(data=data_glomo, aes(x=Sample, y=Abundance, fill=Order))+
  geom_bar(stat="identity", position="stack")+ scale_fill_manual(values = MyOrderPalette) +
  theme_bw() + labs(x = "Sample Name", y = "Relative Abundance (%)") +theme(axis.text.x = element_text(angle = 90))+
  scale_x_discrete(limits=c("3-S58", "37-S92", "2-S57", "4-S59", "25-S80", 
                            "5-S60", "6-S61", "17-S72", "33-S88", 
                            "9-S64", "14-S69", "28-S83", "24-S79", "16-S71", "19-S74", "30-S85", 
                            "1-S56", "15-S70", "20-S75", "23-S78", "22-S77", "21-S76", "26-S81", "31-S86", "27-S82",  "13-S68", "36-S91", 
                            "10-S65", "7-S62", "8-S63", "11-S66", "32-S87", 
                            "12-S67", "18-S73", "29-S84", "34-S89", "35-S90"))
print(OrderGraph)

#Family 1% combo#
MyFamilyPalette <- c("grey", "darkgreen", "steelblue1", "goldenrod1", "violetred1", "cornflowerblue", "dodgerblue", "red", "blue", "chocolate", "peachpuff", "violetred1", "greenyellow", "gold", "violet", "darkslateblue", "cadetblue1", "olivedrab", "palegoldenrod","purple", "yellow3", "plum", "palegreen2", "lightcoral", "chartreuse1", "maroon", "slateblue1", "tomato", "purple", "darkred", "hotpink1", "skyblue1", "darkgreen", "lemonchiffon1", "springgreen1", "maroon2", "dodgerblue", "mediumblue", "yellow3", "coral", "darkolivegreen4", "chartreuse1", "slateblue1", "tomato", "darkmagenta", "purple", "goldenrod1", "lemonchiffon1", "chartreuse1", "maroon", "slateblue1", "tomato", "purple", "darkred", "hotpink1", "skyblue1", "darkgreen", "lemonchiffon1", "springgreen1", "maroon2", "DODGERBLUE", "mediumblue", "yellow3", "coral", "darkolivegreen4", "cornflowerblue", "darkslategrey", "seagreen3", "orange", "darkorchid1", "mistyrose", "cadetblue1", "purple4", "red", "green4", "blue", "plum", "peachpuff", "midnightblue", "indianred3", "khaki", "mediumpurple2", "steelblue1", "greenyellow", "deepskyblue2", "olivedrab", "yellow", "black", "forestgreen", "darkslateblue", "rosybrown1", "grey", "navy", "palegreen")
FamilyGraph <- ggplot(data=data_glomf, aes(x=Sample, y=Abundance, fill=Family))+
  geom_bar(stat="identity", position="stack")+ scale_fill_manual(values = MyFamilyPalette) +
  theme_bw() + labs(x = "Sample Name", y = "Relative Abundance (%)") +theme(axis.text.x = element_text(angle = 90))+
  scale_x_discrete(limits=c("3-S58", "37-S92", "2-S57", "4-S59", "25-S80", 
                            "5-S60", "6-S61", "17-S72", "33-S88", 
                            "9-S64", "14-S69", "28-S83", "24-S79", "16-S71", "19-S74", "30-S85", 
                            "1-S56", "15-S70", "20-S75", "23-S78", "22-S77", "21-S76", "26-S81", "31-S86", "27-S82",  "13-S68", "36-S91", 
                            "10-S65", "7-S62", "8-S63", "11-S66", "32-S87", 
                            "12-S67", "18-S73", "29-S84", "34-S89", "35-S90"))
print(FamilyGraph)


#Genus <2% combo# - mistyrose1 (removed after steelblue1), sienna1 removed before navy
MyGenusPalette <- c("grey", "darkgreen", "steelblue1", "purple", "yellow3", "plum", "palegreen2", "lightcoral", "chartreuse1", "goldenrod1", "violetred1", "cornflowerblue", "dodgerblue", "red", "blue", "chocolate", "peachpuff", "violetred1", "greenyellow", "gold","maroon", "slateblue1", "tomato", "purple", "darkred", "hotpink1", "skyblue1", "darkgreen", "lemonchiffon1", "springgreen1", "maroon2", "dodgerblue", "mediumblue", "yellow3", "coral", "darkolivegreen4", "cornflowerblue", "darkslategrey", "seagreen3", "orange", "darkorchid1", "mistyrose", "cadetblue1", "purple4", "red", "green4", "blue", "plum", "peachpuff", "midnightblue", "indianred3", "khaki", "mediumpurple2", "steelblue1", "greenyellow", "deepskyblue2", "firebrick", "yellow", "lightcoral", "forestgreen", "darkslateblue", "rosybrown1", "grey", "navy", "palegreen", "mistyrose1", "slateblue1", "darkorchid1", "gold", "cadetblue1", "springgreen1", "maroon2", "dodgerblue", "rosybrown1", "darkolivegreen4", "skyblue1", "lavender", "orange","hotpink1", "GREENYELLOW", "LEMONCHIFFON1", "cyan1", "purple4", "red", "green4", "blue", "cornflowerblue", "olivedrab", "mediumpurple1", "darkred", "sandybrown", "navy", "violet", "khaki", "lawngreen", "paleturquoise1", "tan", "forestgreen", "yellow", "black", "tomato", "black", "grey", "sienna1")
GenusGraph <- ggplot(data=data_glomg, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(stat="identity", position="stack")+ 
  scale_fill_manual(values = MyGenusPalette) + 
  theme_bw() + labs(x = "Sample Name", y = "Relative Abundance (%)") +theme(axis.text.x = element_text(angle = 90))+
  scale_x_discrete(limits=c("3-S58", "37-S92", "2-S57", "4-S59", "25-S80", 
                            "5-S60", "6-S61", "17-S72", "33-S88", 
                            "9-S64", "14-S69", "28-S83", "24-S79", "16-S71", "19-S74", "30-S85", 
                            "1-S56", "15-S70", "20-S75", "23-S78", "22-S77", "21-S76", "26-S81", "31-S86", "27-S82",  "13-S68", "36-S91", 
                            "10-S65", "7-S62", "8-S63", "11-S66", "32-S87", 
                            "12-S67", "18-S73", "29-S84", "34-S89", "35-S90"))
print(GenusGraph)




#"3-S58", "37-S92", "2-S57", "4-S59", "25-S80", #ELLESMERE ISLAND
#"5-S60", "6-S61", #EUREKA
#"17-S72", "33-S88", #eureka site 1
#"9-S64", "14-S69", "28-S83", "24-S79", "16-S71", "19-S74", "30-S85", #MOUNT FLEMING
#"1-S56", "15-S70", "20-S75", "23-S78", "22-S77", "21-S76", "26-S81", "31-S86", "27-S82",  "13-S68", "36-S91", #LINNAEUS TERRACE
#"10-S65", "7-S62", "8-S63", #TYROL VALLEY
#"11-S66", "32-S87", #HORSESHOE MOUNTAIN
#"12-S67", "18-S73", "29-S84", "34-S89", "35-S90") #BATTLESHIP



#### Making NMDS Plots #### 

## Merging decontaminated phyloseq file with a mapping file ##
#### PHYLOSEQ FILE FOR CONTINUED USE: ps2.TFS.CT OR ps.WO.C15 ####

library("plyr"); packageVersion("plyr")
theme_set(theme_bw())
#just otu's
k.ord <- ordinate(glomk, "NMDS", "bray")
p1 = plot_ordination(glomk, k.ord, type="taxa", color="Phylum", title="taxa")
print(p1)

#facet wrap
p1 + facet_wrap(~Phylum, 3)
#just samples
p2 = plot_ordination(GP1, GP.ord, type="samples", color="SampleType", shape="human") 
p2 + geom_polygon(aes(fill=SampleType)) + geom_point(size=5) + ggtitle("samples")

#biplot graphic
p3 = plot_ordination(GP1, GP.ord, type="biplot", color="SampleType", shape="Phylum", title="biplot")
# Some stuff to modify the automatic shape scale
GP1.shape.names = get_taxa_unique(GP1, "Phylum")
GP1.shape <- 15:(15 + length(GP1.shape.names) - 1)
names(GP1.shape) <- GP1.shape.names
GP1.shape["samples"] <- 16
p3 + scale_shape_manual(values=GP1.shape)

#split graphic
p4 = plot_ordination(GP1, GP.ord, type="split", color="Phylum", shape="human", label="SampleType", title="split") 
p4


