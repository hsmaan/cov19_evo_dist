library(Biostrings)
library(stringr)
library(DECIPHER)
library(ape)
library(data.table)
library(uwot)

# Functions from CovidGenotyper/R/global.R  

align_get <- function(fastas, align_iters, align_refs, ncores) {
  
  fastas_ungapped <- RemoveGaps(fastas, removeGaps = "all", processors = ncores)
  fasta_align <- AlignSeqs(fastas_ungapped, iterations = align_iters, refinements = align_refs, processors = ncores)
  fasta_mat <- as.matrix(fasta_align)
  fasta_bin <- as.DNAbin(fasta_mat)
  fasta_ungapped <- del.colgapsonly(fasta_bin, threshold = 0.95)
  fasta_string <- fasta_ungapped %>% as.list %>% as.character %>% lapply(., paste0, collapse = "") %>% unlist %>% DNAStringSet
  align_mat <- as.matrix(fasta_string)
  align_mat_sub <- align_mat[, -mask_sites]
  align_mat_bin <- as.DNAbin(align_mat_sub)
  align_masked <- align_mat_bin %>% as.list %>% as.character %>% lapply(., paste0, collapse = "") %>% unlist %>% DNAStringSet
  return(align_masked)
  
}

dist_get <- function(align, metric) {
  
  mask_sites <- c(187, 1059, 2094, 3037, 3130, 6990, 8022, 10323, 10741, 11074, 13408, 14786, 19684, 20148, 21137, 24034, 24378, 25563, 26144, 26461, 26681, 28077, 28826, 28854, 29700, 4050, 13402, 11083, 15324, 21575)
  align_mat <- as.matrix(align)
  align_mat_sub <- align_mat[, -mask_sites]
  align_mat_bin <- as.DNAbin(align_mat_sub)
  align_masked <- align_mat_bin %>% as.list %>% as.character %>% lapply(., paste0, collapse = "") %>% unlist %>% DNAStringSet
  align_trim <- subseq(align_masked, start = 265, end = 29674)
  dec_dist <- dist.dna(as.DNAbin(align_trim), model = metric, as.matrix = TRUE, pairwise.deletion = FALSE)
  colnames(dec_dist) <- (str_split_fixed(colnames(dec_dist), fixed("."), 2)[,1])
  rownames(dec_dist) <- (str_split_fixed(rownames(dec_dist), fixed("."), 2)[,1])
  return(dec_dist)
  
}

# Read args

args <- commandArgs(trailingOnly = TRUE)

fasta_file <- as.character(args[1])
length_cutoff <- as.numeric(args[2])
ambg_cutoff <- as.numeric(args[3])
dist_method <- as.character(args[4])
iters <- as.numeric(args[5])
refs <- as.numeric(args[6])
cores <- as.numeric(args[7])

# Load input fasta

con_fasta <- readDNAStringSet(fasta_file)

# Test if any fastas are <29000 nucleotides in length

len_rm <- (which(width(con_fasta) < length_cutoff))

if (length(len_rm) > 0) {
	con_fasta <- con_fasta[-len_rm]
}

print(paste0(as.character(length(len_rm)), " sequences removed due to nucleotide length < ", as.character(length_cutoff)))

# Remove based on ambg nucleotide freq

ambg_freq_rm <- which(((apply((letterFrequency(con_fasta, c("N", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V"))), 1, sum))/29000) > ambg_cutoff)

if (length(ambg_freq_rm) > 0) {
	con_fasta <- con_fasta[-ambg_freq_rm]
}

print(paste0(as.character(length(ambg_freq_rm)), " sequences removed due to ambiguous nucleotide frequency > ", as.character(ambg_cutoff*100), "%"))

# Get fasta names

fasta_names <- as.list(names(con_fasta))

# Get DNA alignment 

fasta_aligned <- align_get(con_fasta, align_iters = iters, align_refs = refs, ncores = cores)

# Get DNA distance 

fasta_dist <- dist_get(fasta_aligned, metric = dist_method)

