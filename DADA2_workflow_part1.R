#Last updated October 13, 2020
#DADA2 pipeline part 1

library(stringr)
library(dada2)
library(seqRFLP)

#directory to store trimmed reads 
dir.create("./trimmed")
#directory to store filtered reads
dir.create("./filtered")

args <- commandArgs(trailingOnly=TRUE)
reads <- read.table(args[1], sep="\t", encoding = "UTF-8")

both_primers = FALSE #set to TRUE if both primers need to be removed from the fwd and the rev reads (will be the case if the target region is shorter than the read length)

fprimer = args[2]
rprimer = args[3]

min_length = args[4] 
max_length = args[5] 

freads_trimmed <- c()
rreads_trimmed <- c()

freads_filtered <- c()
rreads_filtered <- c()

#remove the primers with cutadapt 
for (row in 1:nrow(reads)){

	sample = reads[row,1]
	fread = reads[row,2]
	rread =  reads[row,3]

	fread_out = paste("./trimmed/", str_split(fread, ".fq", simplify=TRUE)[1], "_trimmed.fastq", sep="")
	freads_trimmed <- c(freads_trimmed, fread_out)
	freads_filtered <- c(freads_filtered, paste("./filtered/", str_split(fread, ".fq", simplify=TRUE)[1], "_trimmed_filtered.fastq", sep=""))

	rread_out = paste("./trimmed/", str_split(rread, ".fq", simplify=TRUE)[1], "_trimmed.fastq", sep="")
	rreads_trimmed <- c(rreads_trimmed, rread_out)
	rreads_filtered <- c(rreads_filtered, paste("./filtered/", str_split(rread, ".fq", simplify=TRUE)[1], "_trimmed_filtered.fastq", sep=""))

	#run cutadapt
	print(paste("Now running cutadapt on sample:", sample))
	if (both_primers){
		#format primer sequences for cutadapt 
		fread_primers = paste(fprimer, ";required...", revComp(rprimer), ";optional", sep="")
		rread_primers = paste(rprimer, ";required...", revComp(fprimer), ";optional", sep="")
		system(paste('cutadapt -a', paste('"', fread_primers, '"', sep=""), '-A', paste('"', rread_primers, '"', sep=""), '-m', min_length, '-M', max_length, '--discard-untrimmed -j 12 -o', fread_out, '-p', rread_out, fread, rread, '>> cutadapt_stats.txt 2>&1'))
	}
	else{
		system(paste('cutadapt -g', fprimer, '-G', rprimer, '-m', min_length, '-M', max_length, '--discard-untrimmed -j 12 -o', fread_out, '-p', rread_out, fread, rread, '>> cutadapt_stats.txt 2>&1'))
	}
}


print("Cutadapt summary:")
system("python /home/kfgibson/scripts/sum_cutadapt.py")

#plot read quality profiles 
pdf("quality_plots.pdf")
plotQualityProfile(freads_trimmed)
plotQualityProfile(rreads_trimmed)
