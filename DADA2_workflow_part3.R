#Last updated October 13, 2020
#DADA2 pipeline part 3

library(stringr)
library(dada2)

set.seed(100) #for reproducibility

seqtab.nochim.pool <- readRDS("seqtab.nochim.pool.rds")

#assign taxonomy to the genus level
#silva
taxa_silva <- assignTaxonomy(seqtab.nochim.pool, "~/scratch/Databases/18s/silva/SILVA_archive_downloads/SILVA_132_SSURef_Nr99_tax_silva_seven_levels_no_species.fasta", multithread=TRUE, minBoot=80)
#pr2
taxa_pr2 <- assignTaxonomy(seqtab.nochim.pool, "~/scratch/Databases/18s/pr2/accessions/pr2_version_4.12.0_18S_dada2_accessions_added_no_species.fa", multithread=TRUE, minBoot=80, taxLevels = c("Kingdom", "Supergroup","Division", "Class", "Order", "Family", "Genus"))


#assign species using exact matches and allow multiple matches
species_silva <- assignSpecies(seqtab.nochim.pool, "~/scratch/Databases/18s/silva/SILVA_archive_downloads/SILVA_132_SSURef_Nr99_tax_silva_seven_levels_species_assign.fasta", allowMultiple=TRUE)

species_pr2 <- assignSpecies(seqtab.nochim.pool, "~/scratch/Databases/18s/pr2/accessions/pr2_version_4.12.0_18S_dada2_accessions_added_species_assign.fa", allowMultiple=TRUE)


#make ASV headers - first make an empty vector
asv_headers <- vector(dim(seqtab.nochim.pool)[2], mode="character")
#assign header names to vector
for (i in 1:dim(seqtab.nochim.pool)[2]){
	asv_headers[i] <- paste(">ASV", i, sep="_")
}

#output #3: taxonomy tables
asv_taxa_silva <- taxa_silva
row.names(asv_taxa_silva) <- sub(">", "", asv_headers)
write.table(asv_taxa_silva, "ASVs_taxonomy_silva.tsv", sep="\t", quote=F, col.names=NA)

asv_taxa_pr2 <- taxa_pr2
row.names(asv_taxa_pr2) <- sub(">", "", asv_headers)
write.table(asv_taxa_pr2, "ASVs_taxonomy_pr2.tsv", sep="\t", quote=F, col.names=NA)

asv_species_silva <- species_silva
row.names(asv_species_silva) <- sub(">", "", asv_headers)
write.table(asv_species_silva, "ASVs_species_silva.tsv", sep="\t", quote=F, col.names=NA)

asv_species_pr2 <- species_pr2
row.names(asv_species_pr2) <- sub(">", "", asv_headers)
write.table(asv_species_pr2, "ASVs_species_pr2.tsv", sep="\t", quote=F, col.names=NA)
