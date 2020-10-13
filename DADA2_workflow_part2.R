#Last updated October 13, 2020
#DADA2 pipeline part 2

library(stringr)
library(dada2)

set.seed(100) #for reproducibility 

getN <- function(x) sum(getUniques(x))

args <- commandArgs(trailingOnly=TRUE)
reads <- read.table(args[1], sep="\t", encoding = "UTF-8")
min_len <- as.numeric(args[2])
trunc_min <-as.numeric(args[3])
trunc_max <- as.numeric(args[4])
min_overlap <- as.numeric(args[5])

freads_trimmed <- c()
rreads_trimmed <- c()

freads_filtered <- c()
rreads_filtered <- c()

samples <- c()

#read in files names
for (row in 1:nrow(reads)){
	
	sample = as.character(reads[row,1])
	samples <- c(samples, sample)
	fread = reads[row,2]
	rread =  reads[row,3]
	
	fread_out = paste("./trimmed/", str_split(fread, ".fq", simplify=TRUE)[1], "_trimmed.fastq", sep="")
	freads_trimmed <- c(freads_trimmed, fread_out)
	freads_filtered <- c(freads_filtered, paste("./filtered/", str_split(fread, ".fq", simplify=TRUE)[1], "_trimmed_filtered.fastq", sep=""))
	
	rread_out = paste("./trimmed/", str_split(rread, ".fq", simplify=TRUE)[1], "_trimmed.fastq", sep="")
	rreads_trimmed <- c(rreads_trimmed, rread_out)
	rreads_filtered <- c(rreads_filtered, paste("./filtered/", str_split(rread, ".fq", simplify=TRUE)[1], "_trimmed_filtered.fastq", sep=""))

}

filtered_out <- filterAndTrim(freads_trimmed, freads_filtered, rreads_trimmed, rreads_filtered, maxEE=c(2,2), rm.phix=TRUE, minLen=min_len, truncLen=c(trunc_min,trunc_max), multithread=12)

#summary of the number of reads in and the number of reads out (after filtration) for each sample$
print("Filtration summary:")
filtered_out

#sanity check; look at the new quality profile plots for the filtered reads$
pdf("filtered_quality_plots.pdf")
plotQualityProfile(freads_filtered)
plotQualityProfile(rreads_filtered)

#learn the errors
err_freads <- learnErrors(freads_filtered, multithread=TRUE)
err_rreads <- learnErrors(rreads_filtered, multithread=TRUE)

#sanity check: look at the graphs for the estimated error rates
pdf("error_plots.pdf")
plotErrors(err_freads, nominalQ=TRUE)
plotErrors(err_rreads, nominalQ=TRUE)

#dereplicate the sequences
derep_fwd <- derepFastq(freads_filtered)
names(derep_fwd) <- samples 
derep_rev <- derepFastq(rreads_filtered)
names(derep_rev) <- samples

#infer ASVs
dada_fwd.pool <- dada(derep_fwd, err=err_freads, multithread=TRUE, pool=TRUE)
dada_rev.pool <- dada(derep_rev, err=err_rreads, multithread=TRUE, pool=TRUE)

#merge reads
merged_amplicons.pool <- mergePairs(dada_fwd.pool, derep_fwd, dada_rev.pool, derep_rev, minOverlap=min_overlap)

#construct ASV count table
seqtab.pool <- makeSequenceTable(merged_amplicons.pool)

#check the size distribtution of the merged sequences 
print("Size distribution of merged sequences:")
table(nchar(getSequences(seqtab.pool)))

#remover chimeras
seqtab.nochim.pool <- removeBimeraDenovo(seqtab.pool, multithread=TRUE)

summary_tab.pool <- data.frame(row.names=samples, input=filtered_out[,1], filtered=filtered_out[,2], dada_f=sapply(dada_fwd.pool, getN), dada_r=sapply(dada_rev.pool, getN), merged=sapply(merged_amplicons.pool, getN), nonchim=rowSums(seqtab.nochim.pool), final_perc_reads_retained=round(rowSums(seqtab.nochim.pool)/filtered_out[,1]*100, 1))
print("Summary:")
summary_tab.pool

#save object to use in taxonomy assignment step
saveRDS(seqtab.nochim.pool, "seqtab.nochim.pool.rds")

#extract ASV sequences
asv_seqs <- colnames(seqtab.nochim.pool)
#make ASV headers - first make an empty vector
asv_headers <- vector(dim(seqtab.nochim.pool)[2], mode="character")
#assign header names to vector
for (i in 1:dim(seqtab.nochim.pool)[2]){
	asv_headers[i] <- paste(">ASV", i, sep="_")
}

#output #1: ASV fasta file 
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

#output #2: count table
asv_tab <- t(seqtab.nochim.pool)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)
