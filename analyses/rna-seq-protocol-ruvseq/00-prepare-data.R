# Merge PBTA, GTEx, TARGET and TCGA RSEM expected read counts.
# 1. Matched `gtex_target_tcga-gene-counts-rsem-expected_count-collapsed.rds`
#    column names with the mapping files shared by
#    [@komalsrathi](https://github.com/komalsrathi) at
#    <https://github.com/PediatricOpenTargets/ticket-tracker/
#         issues/22#issuecomment-854901528>
# 2. Removed duplicated samples. Some TCGA `sample_barcode`s are mapped to
#    multiple `sample_id`s, and some of their expected count sums are different
#    in `gtex_target_tcga-gene-counts-rsem-expected_count-collapsed.rds`.
#    Keep one of the duplicated samples that have the same RESM expected count.
# 3. Use PBTA `gene-counts-rsem-expected_count-collapsed.rds` column names as
#    both `sample_barcode` and `sample_id` for consitency.
# 4. Find the histologies of the samples by matching their `sample_barcode`s
#    with the `Kids_First_Biospecimen_ID`s in `histologies.tsv`.



# load histology df and count df ------------------------------------------
htl_df <- read.delim('../../data/histologies.tsv',
                     stringsAsFactors = FALSE, sep = '\t',
                     header=TRUE)

pb_kf_cnt_df <- readRDS(
    '../../data/gene-counts-rsem-expected_count-collapsed.rds')

gt_ta_tc_cnt_df <- readRDS(
    '../../data/gtex_target_tcga-gene-counts-rsem-expected_count-collapsed.rds')

# input sample ID mapping shared by @komalsrathi at 
# <https://github.com/PediatricOpenTargets/ticket-tracker/issues/22#
#      issuecomment-854900493>
kfbid_colid_mapping_df_list <- list(
    gtex = read.table('input/gtex_mapping.txt', sep = '\t', header = TRUE,
                      stringsAsFactors = FALSE),
    target = read.table('input/target_mapping.txt', sep = '\t', header = TRUE,
                        stringsAsFactors = FALSE),
    tcga = read.table('input/tcga_mapping.txt', sep = '\t', header = TRUE,
                      stringsAsFactors = FALSE)
)
# merge all mapping dfs
kfbid_colid_mapping_mdf <- do.call(rbind, kfbid_colid_mapping_df_list)
rownames(kfbid_colid_mapping_mdf) <- NULL

# Steps:
# 1. select unique sample IDs
# 2. subset expected count matrix
# 3. DGE testing

# merged df has the same number of rows
stopifnot(identical(
    sum(sapply(kfbid_colid_mapping_df_list, nrow)),
    nrow(kfbid_colid_mapping_mdf)))
# sample barcodes are not unique
dup_brcds <- kfbid_colid_mapping_mdf$sample_barcode[
    duplicated(kfbid_colid_mapping_mdf$sample_barcode)]
dup_brcd_df <- kfbid_colid_mapping_mdf[
    kfbid_colid_mapping_mdf$sample_barcode %in% dup_brcds, ]
dup_brcd_df <- dup_brcd_df[order(dup_brcd_df$sample_barcode), ]
dup_brcd_df$rsem_expected_cnt_colSum <- colSums(
    gt_ta_tc_cnt_df[, dup_brcd_df$sample_id])
dup_brcd_df$remove <- duplicated(
    dup_brcd_df[, c('sample_barcode', 'rsem_expected_cnt_colSum')])
dup_rm_sample_ids <- dup_brcd_df$sample_id[dup_brcd_df$remove]

rmdup_kfbid_colid_mapping_mdf <- kfbid_colid_mapping_mdf[
    !(kfbid_colid_mapping_mdf$sample_id %in% dup_rm_sample_ids), ]

# sample_ids are unique
stopifnot(identical(length(kfbid_colid_mapping_mdf$sample_id),
                    length(unique(kfbid_colid_mapping_mdf$sample_id))))

# htl_df Kids_First_Biospecimen_ID are not unique
# duplicated Kids_First_Biospecimen_IDs
dup_bids <- htl_df$Kids_First_Biospecimen_ID[
    duplicated(htl_df$Kids_First_Biospecimen_ID)]
dup_bid_df <- htl_df[htl_df$Kids_First_Biospecimen_ID %in% dup_bids, ]
dup_bid_df <- dup_bid_df[order(dup_bid_df$Kids_First_Biospecimen_ID), ]


# col subset
htl_subset_colnames <- c('Kids_First_Biospecimen_ID', 'experimental_strategy',
                         'sample_type', 'cohort', 'broad_histology',
                         'short_histology', 'cancer_group', 'gtex_group',
                         'gtex_subgroup')

cs_htl_df <- htl_df[, htl_subset_colnames]
# sample_barcode matches the mapping files
colnames(cs_htl_df)[1] <- 'sample_barcode'
# remove duplicated sample_barcode
cs_htl_df <- cs_htl_df[!duplicated(cs_htl_df$sample_barcode), ]
stopifnot(identical(
    nrow(cs_htl_df),
    length(unique(htl_df$Kids_First_Biospecimen_ID))))


# add pbta kf sample in mapping df
stopifnot(all(!(colnames(pb_kf_cnt_df) %in%
                    rmdup_kfbid_colid_mapping_mdf$sample_id)))
stopifnot(all(!(colnames(pb_kf_cnt_df) %in%
                    rmdup_kfbid_colid_mapping_mdf$sample_barcode)))
stopifnot(identical(length(colnames(pb_kf_cnt_df)),
                    length(unique(colnames(pb_kf_cnt_df)))))

pb_kf_colid_df <- data.frame(sample_barcode = colnames(pb_kf_cnt_df),
                             sample_id = colnames(pb_kf_cnt_df),
                             stringsAsFactors = FALSE)

stopifnot(all(!(pb_kf_colid_df$sample_id %in%
                    rmdup_kfbid_colid_mapping_mdf$sample_id)))
stopifnot(all(!(pb_kf_colid_df$sample_barcode %in%
                    rmdup_kfbid_colid_mapping_mdf$sample_barcode)))

m_kfbid_sid_df <- rbind(pb_kf_colid_df, rmdup_kfbid_colid_mapping_mdf)

stopifnot(identical(nrow(m_kfbid_sid_df),
                    length(unique(m_kfbid_sid_df$sample_id))))


m_kfbid_sid_htl_df <- merge(
    m_kfbid_sid_df, cs_htl_df,
    all.x = TRUE, by = 'sample_barcode', sort = TRUE)

stopifnot(identical(sort(m_kfbid_sid_htl_df$sample_id),
                    sort(m_kfbid_sid_df$sample_id)))

# gtex sample ids do not match biospecimen ids


# subset and merge count df -----------------------------------------------
stopifnot(all(colnames(pb_kf_cnt_df) %in% m_kfbid_sid_htl_df$sample_id))
stopifnot(all(colnames(gt_ta_tc_cnt_df) %in%
                  c(dup_rm_sample_ids, m_kfbid_sid_htl_df$sample_id)))

cmn_genes <- intersect(rownames(pb_kf_cnt_df), rownames(gt_ta_tc_cnt_df))

m_cnt_df <- cbind(
    pb_kf_cnt_df[cmn_genes, ],
    gt_ta_tc_cnt_df[
        cmn_genes,
        colnames(gt_ta_tc_cnt_df) %in% m_kfbid_sid_htl_df$sample_id]
)


# save data ---------------------------------------------------------------
# save and load in another R session, in order to reduce memory usage
saveRDS(m_cnt_df,
        '../../scratch/pbta_kf_gtex_target_tcga_rsem_expected_cnt_df.rds')

saveRDS(m_kfbid_sid_htl_df,
        '../../scratch/pbta_kf_gtex_target_tcga_histology_df.rds')
