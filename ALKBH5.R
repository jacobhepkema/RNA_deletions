# ALKBH5 analysis

library(GenomicRanges)
library(Biostrings)
library(Rsamtools)
library(stringr)

source("deletion_utils.R")

options(stringsAsFactors = FALSE)

bam_path <- "output/ALKBH5/SRR9646144_plasmid.bam"
temp_path <- "templates/ALKBH5_mRuby3_plasmid.fa"

templ <- Biostrings::readDNAStringSet(filepath=temp_path)
names(templ) <- "ALKBH5_mRuby3_plasmid"

bam_gr_df <- BAM_to_granges(bam_path)
bam_gr <- makeGRangesFromDataFrame(bam_gr_df)

deletion_df <- get_deletion_df_no_overlaps(bam_gr_df = bam_gr_df, templ = templ, min_length = 14)
deletion_df_gr <- makeGRangesFromDataFrame(deletion_df, keep.extra.columns = TRUE)

genome(deletion_df_gr) <- "ALKBH5_mRuby3_plasmid"
deletion_seqs <- getSeq(templ, deletion_df_gr)
deletion_df$seq <- as.character(deletion_seqs)
deletion_df$seq_diverse <- sapply(deletion_df$seq, FUN=seq_diverse)
deletion_df <- deletion_df[deletion_df$seq_diverse,]

deletion_df$seq[deletion_df$strand == "-"] <- sapply(deletion_df$seq[deletion_df$strand == "-"], rc)
deletion_df$intron_motif <- sapply(deletion_df$seq, intron_motif)

write.csv(deletion_df, "output/ALKBH5/SRR9646144_plasmid_deletions.csv")