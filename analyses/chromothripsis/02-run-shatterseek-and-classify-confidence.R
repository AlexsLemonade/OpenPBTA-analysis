# Laura Egolf (@LauraEgolf) 2021
# Partially adapted from Yang Yang (@yangyangclover) 2020 
  # https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/902e4c49a80dbc69d1c5c9055ddeb5d4f752219d/analyses/sv-analysis/02-shatterseek.R


## ===================== Load Packages =====================
library(ShatterSeek)

# Define Magrittr pipe
`%>%` <- dplyr::`%>%`

# Suppress dplyr summarise output (gets repetitive when looping through samples)
options(dplyr.summarise.inform = FALSE)


## ===================== Define Directory Paths =====================
# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "chromothripsis")
results_dir <- file.path(analysis_dir, "results")
cnv_scratch_dir <- file.path(root_dir, "scratch", "cnv-merged")

# Create results directory if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# Create scratch/cnv-merged/ directory if it doesn't exist
if (!dir.exists(cnv_scratch_dir)) {
  dir.create(cnv_scratch_dir, recursive = TRUE)
}


## ===================== Load Independent Specimen List =====================
independent_specimen_list <- readr::read_tsv(file.path(root_dir, "data", "independent-specimens.wgs.primary-plus.tsv"))

# Create vector with all sample names
bioid <- unique(independent_specimen_list$Kids_First_Biospecimen_ID)


## ===================== Load and Format CNV File =====================
# Read cnv consensus file
cnvconsensus <- readr::read_tsv(file.path(root_dir, "data", "pbta-cnv-consensus.seg.gz"))

# Choose independent specimens 
cnvconsensus <- cnvconsensus %>% 
  dplyr::filter(ID %in% bioid)

# Subset bioid to only samples that have CNV data
  # Note that 20 samples are not included in the CNV consensus because they failed QC 
  # for 2+ callers (see analyses/copy_number_consensus_call/results/uncalled_samples.tsv).
  # So this analysis includes 777 samples instead of the full 797 in the independent
  # specimens list.
bioid <- bioid[bioid %in% cnvconsensus$ID]

# Reformat to fit ShatterSeek input requirements: remove chrY, remove rows with NA copy number, remove "chr" notation
cnvconsensus <- cnvconsensus %>% 
  dplyr::filter(chrom != "chrY", 
                !is.na(cnvconsensus$copy.num)) %>%
  dplyr::mutate(chrom = stringr::str_remove_all(chrom, "chr"))


## ===================== Define function to merge consecutive CN segments =====================

### Explanation:

# ShatterSeek treats every CN segment as a copy number oscillation, even if two consecutive segments 
# share the same CN value. This is true whether or not there are gaps in between consecutive segments 
# (and it can't handle NA values).
# Because the consensus CNV results include gaps, this threw off the ShatterSeek results. Samples with 
# no CN changes (CN=2 for every segment) were called as having copy number oscillations.

# To avoid this problem, this function takes CNV data from a single sample and merges consecutive CN segments 
# that share the same CN value.
# There can be gaps of any length in between the segments - they just need to be consecutive (not adjacent) 
# on the same chromosome and have the same copy number value.

mergeCNsegments <- function(cnv_df) {

  # Order rows by chromosome and start position 
  cnv_df <- dplyr::arrange(cnv_df, chrom, loc.start)
  
  # Define function to increment index_counter if chromosome and CN values for two rows do not match
  compare_adjacent_segs <- function(seg1, seg2) {
    if ((seg1$chrom == seg2$chrom) & 
        (seg1$copy.num == seg2$copy.num)) {
          return(index_counter)
      } else {
          return(index_counter + 1)
      }
    }
  
  # Create index column; loop through rows and increment the index each time the chromosome or CN value 
  # does not match the previous row (starting at second row)
  cnv_df$index <- 0
  index_counter <- 1
  for (row_iter in 2:nrow(cnv_df)) {
    index_counter <- compare_adjacent_segs(cnv_df[row_iter, ], cnv_df[row_iter-1, ])
    cnv_df[row_iter, "index"] <- index_counter
  }
  
  # Separately update index for first row (check whether first row matches second row)
  cnv_df[1, "index"] <- ifelse(compare_adjacent_segs(cnv_df[1,], cnv_df[2,]), 1, 0)
    # Set index to 1 if row 1 matches row 2
    # Set index to 0 if row 1 doesn't match row 2
  
  # Merge rows by selecting minimum loc.start and maximum loc.end for each index value
  # Reorder rows by chromosome and start position
  cnv_df_merged <- cnv_df %>%
    dplyr::group_by(chrom, copy.num, index) %>% 
    dplyr::summarise(loc.start=min(loc.start), loc.end=max(loc.end)) %>% 
    dplyr::arrange(chrom, loc.start)
  
}


## ===================== Run ShatterSeek and Combine Results =====================
count <- 0 # Keep track of sample count (print status while running)
total <- length(bioid) # Total number of samples to run
chromoth_combined <- data.frame() # Data frame to merge summary data from different samples
chromoth_obj_list <- vector("list", total) # List to store ShatterSeek chromoth objects from different samples
names(chromoth_obj_list) <- bioid

# Loop through sample list and run ShatterSeek
for (b in bioid) {
  # Print status
  count=count+1
  print(paste0("Running: ", b, " (", count, " of ", total, ")"))
  
  # Read SV file for current sample
  sv_current <- readr::read_tsv(file.path(root_dir, "scratch", "sv-vcf", paste0(b, "_withoutYandM.tsv")), 
                                col_types = readr::cols(chrom1 = readr::col_character(), # Specify type for problematic columns
                                                        chrom2 = readr::col_character(),
                                                        alt2 = readr::col_character())) 
  
  # Subset CNV dataframe to current sample
  cnv_current <-  cnvconsensus[cnvconsensus$ID == b,]
  
  # If CNV or SV file is empty, jump into next loop
  if (nrow(cnv_current) == 0 | nrow(sv_current) == 0) {
      print(paste0(b," is missing CNV or SV data"))
    next;
  }

  # Merge consecutive CN segments that share the same CN value (see function description)
  cnv_current_merged <- mergeCNsegments(cnv_current)
  
  # Write out merged CNV file (for debugging)
  readr::write_tsv(cnv_current_merged, file.path(cnv_scratch_dir, paste0(b, "_merged_cnv.tsv")))
  
  # Build SV and CNV objects
  SV_data <-
    SVs(
      chrom1 = as.character(sv_current$chrom1),
      pos1 = as.numeric(sv_current$pos1),
      chrom2 = as.character(sv_current$chrom2),
      pos2 = as.numeric(sv_current$pos2),
      SVtype = as.character(sv_current$SVtype),
      strand1 = as.character(sv_current$strand1),
      strand2 = as.character(sv_current$strand2)
    )
  CN_data <-
    CNVsegs(
      chrom = as.character(cnv_current_merged$chrom),
      start = cnv_current_merged$loc.start,
      end = cnv_current_merged$loc.end,
      total_cn = cnv_current_merged$copy.num
    )
  
  # Run ShatterSeek
  chromothripsis <- shatterseek(SV.sample=SV_data,seg.sample=CN_data)

  # Save chromSummary in a combined dataframe
  chromoth_summary <- chromothripsis@chromSummary
  chromoth_summary <- cbind(Kids_First_Biospecimen_ID = b, chromoth_summary) # Add sample ID
  chromoth_combined <- rbind(chromoth_combined, chromoth_summary)
  
  # Save chromoth object in a named list
  chromoth_obj_list[[b]] <- chromothripsis
}

# Notes about ShatterSeek output:
  # ShatterSeek reports details for each chromosome for each sample, even if that chromosome doesn't have CNVs/SVs.
  # ShatterSeek only reports one chromothripsis region per chromosome.


## ===================== Classify high and low confidence chromothripsis regions =====================
# Cutoffs used here are described briefly in ShatterSeek tutorial and in more detail 
# in Cortes-Ciriano et al. (Supplemental Note).
# There is one set of criteria for low confidence chromothripsis regions, and two different 
# sets of criteria for high confidence regions.

### Add FDR correction
chromoth_combined$fdr_fragment_joins <- p.adjust(chromoth_combined$pval_fragment_joins, method = "fdr")
chromoth_combined$fdr_chr_breakpoint_enrichment <- p.adjust(chromoth_combined$chr_breakpoint_enrichment, method = "fdr")
chromoth_combined$fdr_exp_cluster <- p.adjust(chromoth_combined$pval_exp_cluster, method = "fdr")

### Define logic vector for each cutoff

## Low Confidence Cutoff
LC_cutoff <- 
  # At least 6 interleaved intrachromosomal SVs:
  ((chromoth_combined$clusterSize_including_TRA - chromoth_combined$number_TRA) >= 6) &
      # Note: this is equivalent to adding DEL, DUP, h2hINV, and t2tINV
  # At least 4 adjacent segments oscillating between 2 CN states
  (chromoth_combined$max_number_oscillating_CN_segments_2_states >= 4) & 
  # Not significant for the fragment joins test (even distribution of SV types)
  (chromoth_combined$fdr_fragment_joins > 0.2) & 
  # Significant for either the chromosomal enrichment or the exponential distribution of breakpoints test:
  (chromoth_combined$fdr_chr_breakpoint_enrichment < 0.2 | chromoth_combined$fdr_exp_cluster < 0.2)
      # Note: Use pval_exp_cluster, not pval_exp_chr

## High Confidence Cutoff 1
## Note: Same as Low Confidence Cutoff, but more stringent for oscillating CN states
HC_cutoff1 <- 
  # At least 6 interleaved intrachromosomal SVs:
  ((chromoth_combined$clusterSize_including_TRA - chromoth_combined$number_TRA) >= 6) &
  # At least 7 adjacent segments oscillating between 2 CN states:
  (chromoth_combined$max_number_oscillating_CN_segments_2_states >= 7) & 
  # Not significant for the fragment joins test (even distribution of SV types)
  (chromoth_combined$fdr_fragment_joins > 0.2) & 
  # Significant for either the chromosomal enrichment or the exponential distribution of breakpoints test:
  (chromoth_combined$fdr_chr_breakpoint_enrichment < 0.2 | chromoth_combined$fdr_exp_cluster < 0.2)

## High Confidence Cutoff 2
HC_cutoff2 <-
  # At least 3 interleaved intrachromosomal SVs and at least 4 interchromosomal SVs:
  ((chromoth_combined$clusterSize_including_TRA - chromoth_combined$number_TRA) >= 3 & 
     chromoth_combined$number_TRA >= 4) &
  # At least 7 adjacent segments oscillating between 2 CN states:
  (chromoth_combined$max_number_oscillating_CN_segments_2_states >= 7) & 
  # Not significant for the fragment joins test (even distribution of SV types)
  (chromoth_combined$fdr_fragment_joins > 0.2)  

### Annotate each row of ShatterSeek results dataframe with chromothripsis call based on cutoffs for 
### high confidence, low confidence, or either confidence level
# Note "low_conf" calls surpass the low-confidence threshold but *not* the high-confidence threshold
chromoth_combined$call_any_conf <- 0
chromoth_combined[which(LC_cutoff | HC_cutoff1 | HC_cutoff2), "call_any_conf"] <- 1
chromoth_combined$call_high_conf <- 0
chromoth_combined[which(HC_cutoff1 | HC_cutoff2), "call_high_conf"] <- 1
chromoth_combined$call_low_conf <- 0
chromoth_combined[which(LC_cutoff & !(HC_cutoff1 | HC_cutoff2)), "call_low_conf"] <- 1

### Create new dataframe with sample-level summary

# Count the number of chromothripsis regions per sample (high, low, or any confidence)
chromoth_per_sample <- chromoth_combined %>% 
  dplyr::group_by(Kids_First_Biospecimen_ID) %>%
  dplyr::summarize_at(c("call_any_conf", "call_high_conf", "call_low_conf"), sum)
names(chromoth_per_sample) <- c("Kids_First_Biospecimen_ID", 
                                "count_regions_any_conf", "count_regions_high_conf", "count_regions_low_conf")

# Create a summary variable indicating whether a sample has no calls, >=1 low-confidence call, 
# or >=1 high-confidence call (if a sample has both low- and high-confidence calls it will be 
# grouped with high-confidence)
chromoth_per_sample <- chromoth_per_sample %>%
  dplyr::mutate(
    any_regions = dplyr::case_when(
      count_regions_high_conf > 0 ~ "High Confidence",
      count_regions_low_conf > 0 ~ "Low Confidence",
      TRUE ~ "No Calls"),  # Set all others to "No Calls"
    any_regions_logical = (count_regions_any_conf>0)
    )


## =====================  Write out results =====================

# Remove newline characters from column 'inter_other_chroms_coords_all'
# (I think they were inserted for plotting purposes, but they make it hard to read & write the data)
chromoth_combined$inter_other_chroms_coords_all <- gsub("\n", ";", chromoth_combined$inter_other_chroms_coords_all)

# Write out ShatterSeek results with chromothripsis calls defined above
write.table(chromoth_combined, file.path(results_dir, "shatterseek_results_per_chromosome.txt"), 
            sep="\t", quote=F, row.names=F)

# Write out per-sample summary of chromothripsis calls
write.table(chromoth_per_sample, file.path(results_dir, "chromothripsis_summary_per_sample.txt"), 
            sep="\t", quote=F, row.names=F)

# Save list of chromoth objects in scratch directory (used for plotting in script 05)
saveRDS(chromoth_obj_list, file = file.path(root_dir, "scratch", "chromoth_obj_list.rds"))
