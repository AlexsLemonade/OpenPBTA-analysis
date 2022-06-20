# Methods for maf creation for the current release as described in https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/990#issuecomment-819051562
#
# - Normalize the vcf using the above described method
# - Add standard FORMAT fields: ONLY FOR STRELKA2 (GT, AF, AD). This helps for downstream compatibility for some tools
# - VEP annotate
# - Annotate using bcftools to add AF from gnomad (the reference that is used for mutect2)
# - Add soft filter, recommend AF < 0.001 and norm sample DP <= 7
# - Annotate hotspots (protected vcf final)
# - Run vcf2maf (protected maf final)
# - Hard filter selecting PASS or HotSpotAllele=1 (public vcf final)
# - Run vcf2maf on that output for public version (public maf final)
#
#
#
# maf_coltypes is a named vector with column types information for each column in a maf file
# More description about maf format is available
# [here](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/)

library(readr)

maf_coltypes <- cols(Hugo_Symbol=col_character(),
                     Entrez_Gene_Id=col_character(),
                     Center=col_character(),
                     NCBI_Build=col_character(),
                     Chromosome=col_character(),
                     Start_Position=col_double(),
                     End_Position=col_double(),
                     Strand=col_character(),
                     Variant_Classification=col_character(),
                     Variant_Type=col_character(),
                     Reference_Allele=col_character(),
                     Tumor_Seq_Allele1=col_character(),
                     Tumor_Seq_Allele2=col_character(),
                     dbSNP_RS=col_character(),
                     dbSNP_Val_Status=col_logical(),
                     Tumor_Sample_Barcode=col_character(),
                     Matched_Norm_Sample_Barcode=col_character(),
                     Match_Norm_Seq_Allele1=col_character(),
                     Match_Norm_Seq_Allele2=col_character(),
                     Tumor_Sample_UUID=col_character(),
                     Matched_Norm_Sample_UUID=col_character(),
                     HGVSc=col_character(),
                     HGVSp=col_character(),
                     HGVSp_Short=col_character(),
                     Transcript_ID=col_character(),
                     Exon_Number=col_character(),
                     t_depth =col_double(),
                     t_ref_count=col_double(),
                     t_alt_count=col_double(),
                     n_depth=col_double(),
                     n_ref_count=col_double(),
                     n_alt_count=col_double(),
                     all_effects=col_character(),
                     Allele=col_character(),
                     Gene=col_character(),
                     Feature=col_character(),
                     Feature_type=col_character(),
                     Consequence=col_character(),
                     cDNA_position=col_character(),
                     CDS_position=col_character(),
                     Protein_position=col_character(),
                     Amino_acids=col_character(),
                     Codons=col_character(),
                     Existing_variation=col_character(),
                     ALLELE_NUM=col_double(),
                     DISTANCE=col_double(),
                     STRAND_VEP=col_double(),
                     SYMBOL=col_character(),
                     SYMBOL_SOURCE=col_character(),
                     HGNC_ID=col_character(),
                     BIOTYPE=col_character(),
                     CANONICAL=col_character(),
                     CCDS=col_character(),
                     ENSP=col_character(),
                     SWISSPROT=col_character(),
                     TREMBL=col_character(),
                     UNIPARC=col_character(),
                     RefSeq=col_character(),
                     SIFT=col_character(),
                     PolyPhen=col_character(),
                     EXON=col_character(),
                     INTRON=col_character(),
                     DOMAINS=col_character(),
                     AF=col_double(),
                     AFR_AF=col_double(),
                     AMR_AF=col_double(),
                     EAS_AF=col_double(),
                     EUR_AF=col_double(),
                     SAS_AF=col_double(),
                     AA_AF=col_double(),
                     EA_AF=col_double(),
                     CLIN_SIG=col_character(),
                     SOMATIC=col_character(),
                     PUBMED=col_logical(),
                     MOTIF_NAME=col_logical(),
                     MOTIF_POS=col_logical(),
                     HIGH_INF_POS=col_character(),
                     MOTIF_SCORE_CHANGE=col_double(),
                     IMPACT=col_character(),
                     PICK=col_character(),
                     VARIANT_CLASS=col_character(),
                     TSL=col_character(),
                     HGVS_OFFSET=col_character(),
                     PHENO=col_character(),
                     GENE_PHENO=col_character(),
                     FILTER=col_character(),
                     flanking_bps=col_character(),
                     vcf_id=col_character(),
                     vcf_qual=col_double(),
                     gnomAD_AF=col_character(),
                     gnomAD_AFR_AF=col_character(),
                     gnomAD_AMR_AF=col_character(),
                     gnomAD_ASJ_AF=col_character(),
                     gnomAD_EAS_AF=col_character(),
                     gnomAD_FIN_AF=col_character(),
                     gnomAD_NFE_AF=col_character(),
                     gnomAD_OTH_AF=col_character(),
                     gnomAD_SAS_AF=col_character(),
                     vcf_pos=col_number(),
                     HotSpotAllele=col_number()
                     )

# Save as RDS to read into the following scripts while reading maf files
saveRDS(maf_coltypes,"input/maf_coltypes.RDS")

