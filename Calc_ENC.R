library(Biostrings)
library(NeriMisc)
library(data.table)
library(coRdon)
library(dplyr)

# setwd("/home/neri/work/") 

genome_fasta <- readDNAStringSet("Downloads/RiboV1.4/RiboV1.4_Contigs.fasta")
gene_tsv <- fread("Downloads/Simplified_AllORFsInfo.tsv")

# genome_fasta <- readDNAStringSet("RVMT_Zenodo_V4/RiboV1.6_Contigs.fasta")
# gene_tsv <-  fread("RVMT_Zenodo_V4/Tables/Simplified_AllORFsInfo.tsv")
Gene_fasta <- extract_gene_sequences(genome_fasta, gene_tsv)

valid_names <- intersect(names(genome_fasta), gene_tsv$seqid)
Gene_tsv <- gene_tsv[seqid %in% valid_names]
frw <- which(Gene_tsv$strand == "+")
rev <- which(Gene_tsv$strand == "-")
# Now Gene_tsv has the same number of recoreds as Gene_fasta
Gene_tsv$ENC <- -1 
Gene_tsv$ENC[frw] <-  ENC(codonTable(Gene_fasta[frw]))
Gene_tsv$ENC[rev] <-  ENC(codonTable(reverseComplement(Gene_fasta[rev])))
# gene_tsv$ENC <- ENC(codonTable(Gene_fasta[]))
write.csv(Gene_tsv,"./Downloads/Gene_tsv.tsv",sep = '\t')


extract_gene_sequences <- function(genome_fasta, gene_tsv) {
  valid_names <- intersect(names(genome_fasta), gene_tsv$seqid)
  Gene_tsv <- gene_tsv[seqid %in% valid_names]
  
  # Remove any sequence identifiers in genome_fasta that are not present in gene_tsv
  genome_fasta <- genome_fasta[names(genome_fasta) %in% valid_names]
  
  Gene_fasta <- genome_fasta[Gene_tsv$seqid]
  Gene_fasta <- narrow(Gene_fasta, start = Gene_tsv$start, end = Gene_tsv$end)
  names(Gene_fasta) <- Gene_tsv$ORFID
  return(Gene_fasta)
}

ENC_per_seqid <- function(gene_tsv) {
  gene_tsv %>%
    select(seqid, ENC, Length) %>%  # Select only the seqid, ENC and Length columns
    group_by(seqid) %>%
    summarize(weighted_ENC = weighted.mean(ENC, Length)) %>%  # Calculate weighted mean ENC
    distinct(seqid, .keep_all = TRUE) %>%
    select(seqid, ENC = weighted_ENC)  # Rename and select columns
}

meanENC <- ENC_per_seqid(Gene_tsv)


normalizedENC <- ENC_per_seqid(Gene_tsv)



####

ENC_per_genome_table <- distinct(gene_tsv[,c("ENC","seqid")] %>% group_by(seqid)  %>% summarise(
  mean_ENC = mean(ENC),
  seqid = seqid
  #something smart to noramlize per genome
  ))
RemoveNACols()
RemoveNARows()
