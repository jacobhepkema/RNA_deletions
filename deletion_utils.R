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


# Read BAM file as GenomicRanges object. REALLY memory-inefficient. Remind me to not write things like this in R
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


# Given BAM GenomicRanges dataframe and a template, get the span of the reads across the template. 
# Only really makes sense if you align everything to a single sequence
get_spans <- function(bam_gr_df, templ){
  require(GenomicRanges)
  spans <- c()
  pb <- txtProgressBar(0, nrow(bam_gr_df), style=3)
  for(i in 1:nrow(bam_gr_df)){
    setTxtProgressBar(pb, i)
    cigar_numbers <- as.numeric(strsplit(gsub(bam_gr_df$cigar[i], pattern = "\\.", replacement = ""), "[A-Z]+")[[1]])
    cigar_letters <- strsplit(gsub(bam_gr_df$cigar[i], pattern = "\\.", replacement = ""), "[0-9]+")[[1]][-c(1)]
    # remove insertions
    cigar_numbers <- cigar_numbers[cigar_letters != "I"]
    cigar_letters <- cigar_letters[cigar_letters != "I"]
    # remove starting/trailing nucleotides
    cigar_numbers <- cigar_numbers[cigar_letters != "S"]
    cigar_letters <- cigar_letters[cigar_letters != "S"]
    # remove hard clipping
    cigar_numbers <- cigar_numbers[cigar_letters != "H"]
    cigar_letters <- cigar_letters[cigar_letters != "H"]
    spans <- c(spans, sum(cigar_numbers) / width(templ))
  }
  return(spans)
}


# Get deletion information without using overlaps with exons/regions
get_deletion_df_no_overlaps <- function(bam_gr_df,
                                        templ,
                                        min_length = 14){
  require(stringr)
  deletion_names <- c()
  deletion_starts <- c()
  deletion_ends <- c()
  deletion_strands <- c()
  deletion_sizes <- c()
  deletion_idx <- c()
  deletion_span <- c()
  pb <- txtProgressBar(0, nrow(bam_gr_df), style=3)
  for(i in 1:nrow(bam_gr_df)){
    setTxtProgressBar(pb, i)
    cigar_numbers <- as.numeric(strsplit(gsub(bam_gr_df$cigar[i], pattern = "\\.", replacement = ""), "[A-Z]+")[[1]])
    cigar_letters <- strsplit(gsub(bam_gr_df$cigar[i], pattern = "\\.", replacement = ""), "[0-9]+")[[1]][-c(1)]
    # remove insertions
    cigar_numbers <- cigar_numbers[cigar_letters != "I"]
    cigar_letters <- cigar_letters[cigar_letters != "I"]
    # remove starting/trailing nucleotides
    cigar_numbers <- cigar_numbers[cigar_letters != "S"]
    cigar_letters <- cigar_letters[cigar_letters != "S"]
    # remove hard clipping
    cigar_numbers <- cigar_numbers[cigar_letters != "H"]
    cigar_letters <- cigar_letters[cigar_letters != "H"]
    consider <- 1:(length(cigar_numbers))
    # check if read has deletion of at least length X in exon
    if(sum(cigar_numbers[consider][cigar_letters[consider] %in% c("D", "N")] > min_length) != 0){
      # cycle through deletions
      for(d in 1:sum(cigar_letters[consider] %in% c("D", "N"))){
        if((d != 1) & (d != length(cigar_numbers))){
          # check if current deletion has at least length X
          if(cigar_numbers[consider][cigar_letters[consider] %in% c("D", "N")][d] > min_length){
            # get current deletion index
            curr_del_idx <- consider[cigar_letters[consider] %in% c("D", "N")][d]
            # check if CIGAR string before and after deletion cover exonic parts ("M")
            if("M" %in% cigar_letters[consider[1]:(curr_del_idx-1)] &
               "M" %in% cigar_letters[(curr_del_idx+1):consider[length(consider)]]){
              # In case of intron, should not be last within exon
              if((cigar_letters[curr_del_idx] == "N" & curr_del_idx == consider[length(consider)]) | 
                 (cigar_letters[curr_del_idx] == "N" & curr_del_idx == consider[1])){
                next
              } else {
                deletion_names <- c(deletion_names, as.character(bam_gr_df$chrom[i]))
                deletion_start <- bam_gr_df$start[i] + sum(cigar_numbers[1:(max(1, curr_del_idx-1))])
                deletion_end <- bam_gr_df$start[i] + sum(cigar_numbers[1:(max(1, curr_del_idx))]) - 1
                deletion_starts <- c(deletion_starts, deletion_start)
                deletion_ends <- c(deletion_ends, deletion_end)
                deletion_sizes <- c(deletion_sizes, cigar_numbers[consider][cigar_letters[consider] %in% c("D", "N")][d])
                deletion_strands <- c(deletion_strands, as.character(bam_gr_df$strand[i]))
                deletion_idx <- c(deletion_idx, i)
                deletion_span <- c(deletion_span, sum(cigar_numbers) / width(templ))
              }
            }
          }
        }
      }
    }
  }
  if(length(deletion_idx) > 0){
    deletion_gr_df <- data.frame(deletion_idx=deletion_idx,
                                 chrom=deletion_names, 
                                 start=deletion_starts, 
                                 end=deletion_ends, 
                                 strand=deletion_strands, 
                                 size=deletion_sizes,
                                 span=deletion_span)
    deletion_gr_df$region <- paste0(deletion_gr_df$chrom, ":", deletion_gr_df$start, "-", deletion_gr_df$end)
    deletion_gr_df$rname <- bam_gr_df$qname[deletion_gr_df$deletion_idx]
    deletion_gr_df$times_found <- as.numeric(table(deletion_gr_df$region)[deletion_gr_df$region])
    return(deletion_gr_df)
  } else {
    return(-1)
  }
}


# Get deletion information using overlaps with exons/regions
get_deletion_df <- function(bam_gr_df, bam_overlaps, exons_gr,
                            min_length = 14){
  require(stringr)
  exons_starts <- start(exons_gr)
  exons_ends <- end(exons_gr)
  deletion_names <- c()
  deletion_starts <- c()
  deletion_ends <- c()
  deletion_strands <- c()
  deletion_sizes <- c()
  deletion_genes <- c()
  deletion_exon_nrs <- c()
  deletion_query_idx <- c()
  deletion_subject_idx <- c()
  pb <- txtProgressBar(0, length(bam_overlaps), style=3)
  for(i in 1:length(bam_overlaps)){
    setTxtProgressBar(pb, i)
    curr_query_idx <- queryHits(bam_overlaps)[i]
    curr_subject_idx <- subjectHits(bam_overlaps)[i]
    # check if CIGAR string has changes in exonic region
    cigar_numbers <- as.numeric(strsplit(gsub(bam_gr_df$cigar[curr_query_idx], pattern = "\\.", replacement = ""), "[A-Z]+")[[1]])
    cigar_letters <- strsplit(gsub(bam_gr_df$cigar[curr_query_idx], pattern = "\\.", replacement = ""), "[0-9]+")[[1]][-c(1)]
    # remove insertions
    cigar_numbers <- cigar_numbers[cigar_letters != "I"]
    cigar_letters <- cigar_letters[cigar_letters != "I"]
    # remove starting/trailing nucleotides
    cigar_numbers <- cigar_numbers[cigar_letters != "S"]
    cigar_letters <- cigar_letters[cigar_letters != "S"]
    # figure out places within read ranges relative to read start where exons begin and end
    exon_starts_at_pos <- exons_starts[curr_subject_idx] - bam_gr_df$start[curr_query_idx]
    add <- 0
    for(j in 1:length(cigar_numbers)){
      cigar_num <- cigar_numbers[j]
      if(add + cigar_num > exon_starts_at_pos){
        break
      } else {
        add <- add + cigar_num
      }
    }
    exon_ends_at_pos <- exons_ends[curr_subject_idx] - bam_gr_df$start[curr_query_idx]
    add <- 0
    for(q in 1:length(cigar_numbers)){
      cigar_num <- cigar_numbers[q]
      if(add + cigar_num > exon_ends_at_pos){
        break
      } else {
        add <- add + cigar_num
      }
    }
    if(q < j){
      next
    } else {
      consider <- j:q
      # check if read has deletion of at least length X in exon
      if(sum(cigar_numbers[consider][cigar_letters[consider] %in% c("D", "N")] > min_length) != 0){
        # cycle through deletions
        for(d in 1:sum(cigar_letters[consider] %in% c("D", "N"))){
          # check if current deletion has at least length X
          if(cigar_numbers[consider][cigar_letters[consider] %in% c("D", "N")][d] > min_length){
            # get current deletion index
            curr_del_idx <- consider[cigar_letters[consider] %in% c("D", "N")][d]
            # check if deletion is not first or last of CIGAR string
            if(curr_del_idx != 1 & 
               curr_del_idx != length(cigar_numbers)){
              # check if CIGAR string before and after deletion cover exonic parts ("M")
              if("M" %in% cigar_letters[consider[1]:(curr_del_idx-1)] &
                 "M" %in% cigar_letters[(curr_del_idx+1):consider[length(consider)]]){
                # In case of intron, should not be last within exon
                if((cigar_letters[curr_del_idx] == "N" & curr_del_idx == consider[length(consider)]) | 
                   (cigar_letters[curr_del_idx] == "N" & curr_del_idx == consider[1])){
                  next
                }
                deletion_names <- c(deletion_names, as.character(bam_gr_df$chrom[curr_query_idx]))
                deletion_start <- bam_gr_df$start[curr_query_idx] + sum(cigar_numbers[1:(max(1, curr_del_idx-1))])
                deletion_end <- bam_gr_df$start[curr_query_idx] + sum(cigar_numbers[1:(max(1, curr_del_idx))]) - 1
                deletion_starts <- c(deletion_starts, deletion_start)
                deletion_ends <- c(deletion_ends, deletion_end)
                deletion_sizes <- c(deletion_sizes, cigar_numbers[consider][cigar_letters[consider] %in% c("D", "N")][d])
                deletion_strands <- c(deletion_strands, as.character(bam_gr_df$strand[curr_query_idx]))
                deletion_genes <- c(deletion_genes, exons_gr$symbol[curr_subject_idx])
                deletion_exon_nrs <- c(deletion_exon_nrs, exons_gr$exonnumber[curr_subject_idx])
                deletion_query_idx <- c(deletion_query_idx, curr_query_idx)
                deletion_subject_idx <- c(deletion_subject_idx, curr_subject_idx)
              }
            }
          }
        }
      }
    }
  }
  deletion_gr_df <- data.frame(query_idx=deletion_query_idx,
                               subject_idx=deletion_subject_idx,
                               chrom=deletion_names, 
                               start=deletion_starts, 
                               end=deletion_ends, 
                               strand=deletion_strands, 
                               symbol=deletion_genes,
                               exon_nr=deletion_exon_nrs,
                               size=deletion_sizes)
  deletion_gr_df$region <- paste0(deletion_gr_df$chrom, ":", deletion_gr_df$start, "-", deletion_gr_df$end)
  return(deletion_gr_df)
}


intron_motif <- function(seq){ 
  split_seq <- strsplit(seq, "")[[1]]
  return(paste0(paste0(split_seq[c(1,2)], collapse=""), paste0(split_seq[c((length(split_seq)-1), length(split_seq))], collapse=""), collapse=""))
}