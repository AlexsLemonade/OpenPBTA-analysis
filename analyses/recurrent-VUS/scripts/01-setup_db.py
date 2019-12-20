#!/bin/env python

# 01-setup_db.py
#
# Creates and/or fills a database of variant calls.
# Note: requires pandas to be installed, and expects python3
#
# All arguments are optional; only the included tables will be affected.
#
# Arguments:
#   -d DB_FILE, --db-file DB_FILE
#     Path of the database file to use or create. Defaults to `data.sqlite`.
#   --consensus-file CONSENSUS_FILE
#     Path of the MAF-like data file containing the consensus calls(TSV).
#   --meta-file META_FILE, --hist-file META_FILE
#     Path of the metadata/histology data file(TSV).
#   --overwrite           Overwrite tables that may already exist.
#
#
# Example invocation:
# python3 01-setup_db.py --db-file snv.sqlite --strelka-file strelka.maf --meta-file samples.tsv

import argparse
import base64
import os
import requests
import sqlite3

import pandas as pd

COSMIC_URL = "https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v90/CosmicMutantExport.tsv.gz"

parser = argparse.ArgumentParser(
    description="Creates and/or fills a database of variant calls")

parser.add_argument(
    '-d',
    '--db-file',
    dest='db_file',
    default='data.sqlite',
    help="Path of the database file to use or create. Defaults to `data.sqlite`."
)

parser.add_argument(
    '--consensus-file',
    dest='consensus_file',
    help="Path of the MAF-like TSV file with consensus SNV calls."
)
parser.add_argument(
    '--meta-file',
    '--hist-file',
    dest='meta_file',
    help="Path of the metadata/histology data file (TSV)."
)
parser.add_argument(
    '--ind-file',
    dest='ind_file',
    help="Path of the independent samples file."
)
parser.add_argument(
    '--overwrite',
    action='store_true',
    help="Flag for whether to overwrite tables that may already exist."
)

# parser.add_argument(
#     '--cosmic-user',
#     dest='cosmic_user',
#     help="Username (email) for COSMIC user."
# )
# 
# parser.add_argument(
#     '--cosmic-pass',
#     dest='cosmic_pass',
#     help="password for COSMIC user."
# )


args = parser.parse_args()

# types for all expected MAF fields.
maf_types = [
    ('Hugo_Symbol', 'TEXT'),
    ('Entrez_Gene_Id', 'INTEGER'),
    ('Center', 'TEXT'),
    ('NCBI_Build', 'TEXT'),
    ('Chromosome', 'TEXT'),
    ('Start_Position', 'INTEGER'),
    ('End_Position', 'INTEGER'),
    ('Strand', 'TEXT'),
    ('Variant_Classification', 'TEXT'),
    ('Variant_Type', 'TEXT'),
    ('Reference_Allele', 'TEXT'),
    ('Tumor_Seq_Allele1', 'TEXT'),
    ('Tumor_Seq_Allele2', 'TEXT'),
    ('dbSNP_RS', 'TEXT'),
    ('dbSNP_Val_Status', 'TEXT'),
    ('Tumor_Sample_Barcode', 'TEXT'),
    ('Matched_Norm_Sample_Barcode', 'TEXT'),
    ('Match_Norm_Seq_Allele1', 'TEXT'),
    ('Match_Norm_Seq_Allele2', 'TEXT'),
    ('Tumor_Validation_Allele1', 'TEXT'),
    ('Tumor_Validation_Allele2', 'TEXT'),
    ('Match_Norm_Validation_Allele1', 'TEXT'),
    ('Match_Norm_Validation_Allele2', 'TEXT'),
    ('Verification_Status', 'TEXT'),
    ('Validation_Status', 'TEXT'),
    ('Mutation_Status', 'TEXT'),
    ('Sequencing_Phase', 'TEXT'),
    ('Sequence_Source', 'TEXT'),
    ('Validation_Method', 'TEXT'),
    ('Score', 'TEXT'),
    ('BAM_File', 'TEXT'),
    ('Sequencer', 'TEXT'),
    ('Tumor_Sample_UUID', 'TEXT'),
    ('Matched_Norm_Sample_UUID', 'TEXT'),
    ('HGVSc', 'TEXT'),
    ('HGVSp', 'TEXT'),
    ('HGVSp_Short', 'TEXT'),
    ('Transcript_ID', 'TEXT'),
    ('Exon_Number', 'TEXT'),
    ('t_depth', 'INTEGER'),
    ('t_ref_count', 'INTEGER'),
    ('t_alt_count', 'INTEGER'),
    ('n_depth', 'INTEGER'),
    ('n_ref_count', 'INTEGER'),
    ('n_alt_count', 'INTEGER'),
    ('all_effects', 'TEXT'),
    ('Allele', 'TEXT'),
    ('Gene', 'TEXT'),
    ('Feature', 'TEXT'),
    ('Feature_type', 'TEXT'),
    ('Consequence', 'TEXT'),
    ('cDNA_position', 'TEXT'),
    ('CDS_position', 'TEXT'),
    ('Protein_position', 'TEXT'),
    ('Amino_acids', 'TEXT'),
    ('Codons', 'TEXT'),
    ('Existing_variation', 'TEXT'),
    ('ALLELE_NUM', 'INTEGER'),
    ('DISTANCE', 'INTEGER'),
    ('STRAND_VEP', 'TEXT'),
    ('SYMBOL', 'TEXT'),
    ('SYMBOL_SOURCE', 'TEXT'),
    ('HGNC_ID', 'TEXT'),
    ('BIOTYPE', 'TEXT'),
    ('CANONICAL', 'TEXT'),
    ('CCDS', 'TEXT'),
    ('ENSP', 'TEXT'),
    ('SWISSPROT', 'TEXT'),
    ('TREMBL', 'TEXT'),
    ('UNIPARC', 'TEXT'),
    ('RefSeq', 'TEXT'),
    ('SIFT', 'TEXT'),
    ('PolyPhen', 'TEXT'),
    ('EXON', 'TEXT'),
    ('INTRON', 'TEXT'),
    ('DOMAINS', 'TEXT'),
    ('AF', 'TEXT'),
    ('AFR_AF', 'TEXT'),
    ('AMR_AF', 'TEXT'),
    ('ASN_AF', 'TEXT'),
    ('EAS_AF', 'TEXT'),
    ('EUR_AF', 'TEXT'),
    ('SAS_AF', 'TEXT'),
    ('AA_AF', 'TEXT'),
    ('EA_AF', 'TEXT'),
    ('CLIN_SIG', 'TEXT'),
    ('SOMATIC', 'TEXT'),
    ('PUBMED', 'TEXT'),
    ('MOTIF_NAME', 'TEXT'),
    ('MOTIF_POS', 'INTEGER'),
    ('HIGH_INF_POS', 'TEXT'),
    ('MOTIF_SCORE_CHANGE', 'REAL'),
    ('IMPACT', 'TEXT'),
    ('PICK', 'TEXT'),
    ('VARIANT_CLASS', 'TEXT'),
    ('TSL', 'TEXT'),
    ('HGVS_OFFSET', 'TEXT'),
    ('PHENO', 'TEXT'),
    ('MINIMISED', 'TEXT'),
    ('ExAC_AF', 'TEXT'),
    ('ExAC_AF_AFR', 'TEXT'),
    ('ExAC_AF_AMR', 'TEXT'),
    ('ExAC_AF_EAS', 'TEXT'),
    ('ExAC_AF_FIN', 'TEXT'),
    ('ExAC_AF_NFE', 'TEXT'),
    ('ExAC_AF_OTH', 'TEXT'),
    ('ExAC_AF_SAS', 'TEXT'),
    ('GENE_PHENO', 'TEXT'),
    ('FILTER', 'TEXT'),
    ('flanking_bps', 'TEXT'),
    ('vcf_id', 'TEXT'),
    ('vcf_qual', 'REAL'),
    ('ExAC_AF_Adj', 'TEXT'),
    ('ExAC_AC_AN_Adj', 'TEXT'),
    ('ExAC_AC_AN', 'TEXT'),
    ('ExAC_AC_AN_AFR', 'TEXT'),
    ('ExAC_AC_AN_AMR', 'TEXT'),
    ('ExAC_AC_AN_EAS', 'TEXT'),
    ('ExAC_AC_AN_FIN', 'TEXT'),
    ('ExAC_AC_AN_NFE', 'TEXT'),
    ('ExAC_AC_AN_OTH', 'TEXT'),
    ('ExAC_AC_AN_SAS', 'TEXT'),
    ('ExAC_FILTER', 'TEXT'),
    ('gnomAD_AF', 'TEXT'),
    ('gnomAD_AFR_AF', 'TEXT'),
    ('gnomAD_AMR_AF', 'TEXT'),
    ('gnomAD_ASJ_AF', 'TEXT'),
    ('gnomAD_EAS_AF', 'TEXT'),
    ('gnomAD_FIN_AF', 'TEXT'),
    ('gnomAD_NFE_AF', 'TEXT'),
    ('gnomAD_OTH_AF', 'TEXT'),
    ('gnomAD_SAS_AF', 'TEXT'),
    ('vcf_pos', 'INTEGER'),
    ('VAF', 'REAL')
]

# indexes to create. lists to handle multiple column indexes.
maf_indexes = [
    ('SNV', ['Chromosome', 'Start_Position', 'Reference_Allele', 'Allele']),
    ('position', ['Chromosome', 'Start_Position']),
    ('sample', ['Tumor_Sample_Barcode'])
]

# translate sqltypes to dtypes
dtypes = {
    'TEXT': 'object',
    'INTEGER': 'Int64',
    'REAL': 'float64'
}
maf_dtypes = {col: dtypes[sql_type] for col, sql_type in maf_types}

# connect to/create the db
con = sqlite3.connect(args.db_file)

chunksize = 1e5

####### Create consensus SNV table
print("Reading file {} to table 'consensus'.".format(args.consensus_file))
# we need 2 full tables for consensus with one missing
if args.overwrite:
    con.execute("DROP TABLE IF EXISTS consensus")
maf_table_def = ", ".join([" ".join(col) for col in maf_types])
maf_create = "CREATE TABLE consensus('index' INTEGER, {})".format(maf_table_def)
con.execute(maf_create)
# read and fill maf table
maf_chunks = pd.read_table(args.consensus_file,
                           dtype=maf_dtypes,
                           na_values=["."],
                           comment="#",
                           chunksize=chunksize)

for chunk in maf_chunks:
    # process the chunk
    chunk['VAF'] = (chunk['t_alt_count'] /
                    (chunk['t_ref_count'] + chunk['t_alt_count']))
    chunk.to_sql('consensus', con, if_exists='append')
# create indexes
print("Indexing table 'consensus'")
for index_name, fields in maf_indexes:
    index_statement = "CREATE INDEX consensus_{name} ON consensus({field})".format(
        name=index_name,
        field=", ".join(fields)
    )
    con.execute(index_statement)

######### Read the metadata file and load into a table called "samples"
if args.meta_file:
    print("Reading file {} to table 'samples'".format(args.meta_file))
    if args.overwrite:
        con.execute("DROP TABLE IF EXISTS samples")
    metadata_df = pd.read_table(args.meta_file)
    metadata_df.to_sql("samples", con)

if args.ind_file:
    print("Reading file {} to table 'ind_samples'".format(args.ind_file))
    if args.overwrite:
        con.execute("DROP TABLE IF EXISTS ind_samples")
    metadata_df = pd.read_table(args.ind_file)
    metadata_df.to_sql("ind_samples", con)
