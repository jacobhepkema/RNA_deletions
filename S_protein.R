# S_protein analysis

suppressMessages(library(stringr))
suppressMessages(library(GenomicRanges))
suppressMessages(library(Biostrings))
suppressMessages(library(BSgenome))
suppressMessages(library(Rsamtools))
suppressMessages(source("deletion_utils.R"))
options(stringsAsFactors = FALSE)

ID_samp_temp <- read.csv("data/S_protein_samples_templates.csv")

for(row_idx in 1:nrow(ID_samp_temp)){
  temp_path <- ID_samp_temp$template[row_idx]
  templ <- Biostrings::readDNAStringSet(filepath=temp_path)
  names(templ) <- str_remove(basename(temp_path), "\\.fa")
  
  ID <- ID_samp_temp$ID[row_idx]
  cat(paste0("\n", ID, "\n"))
  
  bam_path <- paste0("output/S_protein/", ID, ".bam")
  bam_gr_df <- BAM_to_granges(bam_path)
  bam_gr <- makeGRangesFromDataFrame(bam_gr_df)
  
  cat("\n")
  
  deletion_df <- get_deletion_df_no_overlaps(bam_gr_df = bam_gr_df, templ = templ, min_length = 14)
  deletion_df_gr <- makeGRangesFromDataFrame(deletion_df, keep.extra.columns = TRUE)
  
  genome(deletion_df_gr) <- str_remove(basename(temp_path), "\\.fa")
  deletion_seqs <- Biostrings::getSeq(templ, deletion_df_gr)
  deletion_df$seq <- as.character(deletion_seqs)
  deletion_df$seq_diverse <- sapply(deletion_df$seq, FUN=seq_diverse)
  deletion_df <- deletion_df[deletion_df$seq_diverse,]
  
  deletion_df$seq[deletion_df$strand == "-"] <- sapply(deletion_df$seq[deletion_df$strand == "-"], rc)
  deletion_df$intron_motif <- sapply(deletion_df$seq, intron_motif)
  
  write.csv(deletion_df, paste0("output/S_protein/", ID, "_deletions.csv"))
}

print("Done")
