library(data.table) #used for data table functions
library(ggplot2) #used for charts plotting

setwd(".")
# use fread since we have large file
dt_orig_output <- fread(paste0(getwd(),"/mrbwa_output.tsv"), sep = '\t', header = FALSE)
colnames(dt_orig_output) <- c('ref_index','ref_bp','chromosone','read_bps_counts')

read_bps <- c('A','C','G','T','D','N')

# initialise final data frame with colnames
df_analysed_data <- data.frame(ref_index=NA,
                               ref_bp=NA,
                               chromosone=NA,
                               read_bp=NA,
                               read_bp_correctness=NA,
                               matching=NA,
                               insertion=NA,
                               deletion=NA)[numeric(0), ]

# colnames(df_analysed_data) <- c('ref_index','ref_bp','chromosone','read_bp','read_bp_correctness','matching','insertion','deletion')

for (i in 1:nrow(dt_orig_output)) {
  aligned_row <- dt_orig_output[i]
  # parse query counts and extract the aligned bp
  read_bps_counts <- as.numeric(unlist(strsplit(as.character(aligned_row[, 'read_bps_counts', with=FALSE]), ",")))
  read_bps_total_count <- sum(read_bps_counts)
  # create mini data table having bp, count, percentage
  dt_read_bps <- data.table(bp=read_bps, count=seq(0,0, length.out = 6), ratio=seq(0,0, length.out = 6))
  for (j in seq_along(read_bps_counts)) {
    dt_read_bps[j]$count <- read_bps_counts[j]
    dt_read_bps[j]$ratio <- read_bps_counts[j] / read_bps_total_count
  }
  # sort data table by count desc
  dt_read_bps <- dt_read_bps[order(-count)]
  
  read_bp <- dt_read_bps[1]$bp
  read_bp_correctness = dt_read_bps[1]$ratio
  ref_bp <- aligned_row[, 'ref_bp', with=FALSE]
  matching <- ref_bp == read_bp
  
  # we can have 2 bps that have practically the same count and are still correct.
  # This happens because of zygosity i.e. 
  # humans inherit set of chromosones from father and another from mother. hence homologous chromosones can have same or different alleles
  # Hence if read bp does not match, we check if the next bp in line is within a predefined threshold from the previous one and also has at least 30% of the reads
  
  # threshold is defined as the reciprocal of the total reads in order to cater for small and bug number of reads
  zygosity_threshold <- round(1 / read_bps_total_count, 2)
  next_bp_counter <- 2
  next_bp_ratio <- dt_read_bps[next_bp_counter]$ratio
  while (!matching && next_bp_ratio > 0.3 && (next_bp_counter <= nrow(dt_read_bps))) {
    if (read_bp_correctness - next_bp_ratio <= zygosity_threshold) {
      read_bp <- dt_read_bps[next_bp_counter]$bp
      read_bp_correctness <- dt_read_bps[next_bp_counter]$ratio
      
      matching <- ref_bp == read_bp
    }
    if (!matching) {
      next_bp_counter <- next_bp_counter + 1
      # we check for length due to a bug related to comparing with numeric(0) - https://stat.ethz.ch/pipermail/r-help/2008-November/180839.html
      next_bp_ratio <- if (length(dt_read_bps[next_bp_ratio]$ratio) == 0) 0 else dt_read_bps[next_bp_ratio]$ratio
    }
  }
  
  # if still not matching, revert values to original
  if (!matching) {
    read_bp <- dt_read_bps[1]$bp
    read_bp_correctness = dt_read_bps[1]$ratio
  }
  
  df_analysed_data <- rbind(df_analysed_data, 
                            data.frame(aligned_row[, 'ref_index', with=FALSE],
                                       ref_bp,
                                       aligned_row[, 'chromosone', with=FALSE],
                                       read_bp,
                                       read_bp_correctness = round(read_bp_correctness, 2),
                                       matching = ref_bp == read_bp,
                                       insertion = ref_bp == 'None',
                                       deletion = read_bp == 'D'))
}

dt_analysed_data <- data.table(df_analysed_data)

# create data table and plot containing alignment operation counts i.e. match / variation / insertion / deletion
dt_alignment_dist <- data.table(
  Operation=c('Match','Variation','Insertion','Deletion'),
  Count=c(nrow(dt_analysed_data[dt_analysed_data$ref_bp.1 & !dt_analysed_data$deletion & !dt_analysed_data$ref_bp.2]),
          nrow(dt_analysed_data[!dt_analysed_data$ref_bp.1 & !dt_analysed_data$deletion & !dt_analysed_data$ref_bp.2]),
          nrow(dt_analysed_data[dt_analysed_data$ref_bp.2]),
          nrow(dt_analysed_data[dt_analysed_data$deletion])))
# factor for ordering
dt_alignment_dist$Operation=factor(dt_alignment_dist$Operation,levels=dt_alignment_dist$Operation)
# plot
dt_alignment_dist.plot <- ggplot(dt_alignment_dist, aes(x=factor(Operation), y=Count)) +
  geom_bar(stat="identity") + labs(x='Alignment Operation', y='Count') + 
  ggtitle("SRR062634: Alignment Operation Counts") +
  theme(text=element_text(size=15))

# save plot
ggsave(filename=paste0(getwd(), "/alignment_dist.jpg"), plot=dt_alignment_dist.plot)

# mismatch per chromosome. a mismatch is neither a deletion or insertion so we need to filter these out
dt_mismatches <- dt_analysed_data[!dt_analysed_data$ref_bp.1 & !dt_analysed_data$deletion & !dt_analysed_data$ref_bp.2]
dt_mismatches_per_chromosone <- dt_mismatches[, .N, by = chromosone]
colnames(dt_mismatches_per_chromosone) <- c('Chromosone','Count')
# keep only the top mismatches i.e. having > count >= 3rd quartile
dt_mismatches_per_chromosone <- dt_mismatches_per_chromosone[Count >= quantile(dt_mismatches_per_chromosone$Count)[[4]]]
dt_mismatches_per_chromosone <- dt_mismatches_per_chromosone[order(-Count)]
# factor for ordering
dt_mismatches_per_chromosone$Chromosone=factor(dt_mismatches_per_chromosone$Chromosone,levels=dt_mismatches_per_chromosone$Chromosone)
# plot
dt_mismatches_per_chromosone.plot <- ggplot(dt_mismatches_per_chromosone, aes(x=factor(Chromosone), y=Count)) +
  geom_bar(stat="identity") + labs(x='Chromosome', y='Count') + 
  ggtitle("SRR062634: No of Variations per Chromosome") +
  theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust = 1, size = 10))

# save plot
ggsave(filename=paste0(getwd(), "/variations_per_chromosome.jpg"), plot=dt_mismatches_per_chromosone.plot)

# insertions per chromosome
dt_insertions <- dt_analysed_data[dt_analysed_data$ref_bp.2]
dt_insertions_per_chromosone <- dt_insertions[, .N, by = chromosone]
colnames(dt_insertions_per_chromosone) <- c('Chromosone','Count')
# keep only the top insertions i.e. having > count >= 2nd quartile i.e. mean
dt_insertions_per_chromosone <- dt_insertions_per_chromosone[Count >= quantile(dt_insertions_per_chromosone$Count)[[3]]]
dt_insertions_per_chromosone <- dt_insertions_per_chromosone[order(-Count)]
# factor for ordering
dt_insertions_per_chromosone$Chromosone=factor(dt_insertions_per_chromosone$Chromosone,levels=dt_insertions_per_chromosone$Chromosone)
# plot
dt_insertions_per_chromosone.plot <- ggplot(dt_insertions_per_chromosone, aes(x=factor(Chromosone), y=Count)) +
  geom_bar(stat="identity") + labs(x='Chromosome', y='Count') + 
  ggtitle("SRR062634: No of Insertions per Chromosome") +
  theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust = 1, size = 10))

# save plot
ggsave(filename=paste0(getwd(), "/insertions_per_chromosome.jpg"), plot=dt_insertions_per_chromosone.plot)

# deletion per chromosome
dt_deletions <- dt_analysed_data[dt_analysed_data$deletion]
dt_deletions_per_chromosone <- dt_deletions[, .N, by = chromosone]
colnames(dt_deletions_per_chromosone) <- c('Chromosone','Count')
# keep only the top insertions i.e. having > count >= 2nd quartile i.e. mean
dt_deletions_per_chromosone <- dt_deletions_per_chromosone[Count >= quantile(dt_deletions_per_chromosone$Count)[[3]]]
dt_deletions_per_chromosone <- dt_deletions_per_chromosone[order(-Count)]
# factor for ordering
dt_deletions_per_chromosone$Chromosone=factor(dt_deletions_per_chromosone$Chromosone,levels=dt_deletions_per_chromosone$Chromosone)
# plot
dt_deletions_per_chromosone.plot <- ggplot(dt_deletions_per_chromosone, aes(x=factor(Chromosone), y=Count)) +
  geom_bar(stat="identity") + labs(x='Chromosome', y='Count') + 
  ggtitle("SRR062634: No of Deletions per Chromosome") +
  theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust = 1, size = 10))

# save plot
ggsave(filename=paste0(getwd(), "/deletions_per_chromosome.jpg"), plot=dt_deletions_per_chromosone.plot)