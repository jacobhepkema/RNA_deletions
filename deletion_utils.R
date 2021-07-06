# Reverse complement sequences in list
rc <- function(seq_list){
  require(stringr)
  COMPS <- list("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  revcomps <- c()
  for(i in str_split(seq_list, "")){
    revcomps <- c(revcomps, paste0(unlist(COMPS[rev(i)]), collapse=""))
  }
  return(revcomps)
}


# Figure out if fraction of two nucleotides is not over a certain fraction of the sequence
seq_diverse <- function(seq, frac=0.8){
  require(stringr)
  split_seq <- strsplit(seq, "")[[1]]
  a_frac <- sum(split_seq==c("A")) / str_length(seq)
  c_frac <- sum(split_seq==c("C")) / str_length(seq)
  g_frac <- sum(split_seq==c("G")) / str_length(seq)
  t_frac <- sum(split_seq==c("T")) / str_length(seq)
  if(a_frac + c_frac > frac |
     a_frac + g_frac > frac |
     a_frac + t_frac > frac |
     c_frac + g_frac > frac |
     c_frac + t_frac > frac |
     g_frac + t_frac > frac){
    return(FALSE)
  } else {
    return(TRUE)
  }
}


# Read BAM file as GenomicRanges object. REALLY memory-inefficient
BAM_to_granges <- function(path, 
                           param = ScanBamParam(what=c("qname", "rname", "strand", "pos", "qwidth", "cigar", "seq", "mapq")),
                           min_mapq = 10, 
                           remove_chimeric = TRUE){
  require(stringr)
  require(GenomicRanges)
  require(Rsamtools)
  bamfile <- scanBam(path)[[1]]
  good_idx <- which(bamfile$mapq > min_mapq)
  gr_r_names <- bamfile$qname[good_idx]
  gr_names <- bamfile$rname[good_idx]
  gr_start <- bamfile$pos[good_idx]
  gr_strand <- bamfile$strand[good_idx]
  gr_cigar <- bamfile$cigar[good_idx]
  gr_seq <- as.character(bamfile$seq[good_idx])
  bam_gr_df <- data.frame(chrom=gr_names, start=gr_start, end=rep(0,length(gr_start)), # end placeholder
                          strand=gr_strand, cigar=gr_cigar, seq=gr_seq, qname=gr_r_names)
  # I once saw some CIGAR strings with dots in them. Bit strange. Hence the gsub below
  bam_gr_df$mapwidth <- mapply(bam_gr_df$cigar, FUN=function(cigar_string){
    cigar_numbers <- as.numeric(strsplit(gsub(cigar_string, pattern = "\\.", replacement = ""), "[A-Z]+")[[1]])
    cigar_letters <- strsplit(gsub(cigar_string, pattern = "\\.", replacement = ""), "[0-9]+")[[1]][-c(1)]
    # remove insertions
    cigar_numbers <- cigar_numbers[cigar_letters != "I"]
    cigar_letters <- cigar_letters[cigar_letters != "I"]
    # remove trailing nucleotides on both ends
    cigar_numbers <- cigar_numbers[cigar_letters != "S"]
    cigar_letters <- cigar_letters[cigar_letters != "S"]
    return(sum(cigar_numbers))
  })
  bam_gr_df$end <- bam_gr_df$start + bam_gr_df$mapwidth
  if(remove_chimeric){
    old_size <- nrow(bam_gr_df)
    duplicated_read_names <- bam_gr_df$qname[duplicated(bam_gr_df$qname)]
    bam_gr_df <- bam_gr_df[!(bam_gr_df$qname %in% duplicated_read_names),]
    cat(paste0("Removed ", old_size - nrow(bam_gr_df), " chimeric reads\n"))
  }
  return(bam_gr_df)
}
