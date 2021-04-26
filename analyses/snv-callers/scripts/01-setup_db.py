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
#   --strelka-file STRELKA_FILE
#     Path of the MAF formatted data file from the strelka2 caller(TSV).
#   --mutect-file MUTECT_FILE
#     Path of the MAF formatted data file from the mutect2 caller(TSV).
#   --lancet-file LANCET_FILE
#     Path of the MAF formatted data file from the lancet caller(TSV).
#   --vardict-file VARDICT_FILE
#     Path of the MAF formatted data file from the vardict caller(TSV).
#   --meta-file META_FILE, --hist-file META_FILE
#     Path of the metadata/histology data file(TSV).
#   --overwrite           Overwrite tables that may already exist.
#
#
# Example invocation:
# python3 01-setup_db.py --db-file snv.sqlite --strelka-file strelka.maf --meta-file samples.tsv


import os
import sqlite3
import argparse
import pandas as pd

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
    '--strelka-file',
    dest='strelka_file',
    help="Path of the MAF formatted data file from the strelka2 caller (TSV)."
)
parser.add_argument(
    '--mutect-file',
    dest='mutect_file',
    help="Path of the MAF formatted data file from the mutect2 caller (TSV)."
)
parser.add_argument(
    '--lancet-file',
    dest='lancet_file',
    help="Path of the MAF formatted data file from the lancet caller (TSV)."
)
parser.add_argument(
    '--vardict-file',
    dest='vardict_file',
    help="Path of the MAF formatted data file from the vardict caller (TSV)."
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
    '--overwrite',
    action='store_true',
    help="Flag for whether to overwrite tables that may already exist."
)
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
    ('GENE_PHENO', 'TEXT'),
    ('FILTER', 'TEXT'),
    ('flanking_bps', 'TEXT'),
    ('vcf_id', 'TEXT'),
    ('vcf_qual', 'REAL'),
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
    ('HotSpotAllele', 'INTEGER'),
    ('VAF', 'REAL')
]

common_cols = [col for col, type in maf_types]

needed_cols = [
    'Hugo_Symbol',
    'Entrez_Gene_Id',
    'Chromosome',
    'Start_Position',
    'End_Position',
    'Variant_Classification',
    'Variant_Type',
    'Reference_Allele',
    'Tumor_Seq_Allele1',
    'Tumor_Seq_Allele2',
    'dbSNP_RS',
    'Tumor_Sample_Barcode',
    'Transcript_ID',
    't_depth',
    't_ref_count',
    't_alt_count',
    'n_depth',
    'n_ref_count',
    'n_alt_count',
    'Allele',
    'flanking_bps',
    'Protein_position',
    'IMPACT',
    'HGVSp_Short',
    'gnomAD_AF',
    'VAF'
]

needed_types = [col for col in maf_types if col[0] in needed_cols]

# indexes to create. lists to handle multiple column indexes.
maf_indexes = [
    ('SNV', ['Chromosome', 'Start_Position', 'Reference_Allele', 'Allele']),
    ('LOC', ['Chromosome', 'Start_Position']),
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

# Create caller list
callers = [
    ('strelka', args.strelka_file),
    ('mutect', args.mutect_file),
    ('lancet', args.lancet_file),
    ('vardict', args.vardict_file),
    ('consensus', args.consensus_file)
]

for table_name, maf_file in callers:
    if not maf_file:
        continue
    print("Reading file {} to table '{}'.".format(maf_file, table_name, ))
    # we need 2 full tables for consensus with one missing
    if table_name in ('strelka', 'lancet', 'consensus'):
        table_types = maf_types
    else:
        table_types = needed_types
    if args.overwrite:
        con.execute("DROP TABLE IF EXISTS {}".format(table_name))
    maf_table_def = ", ".join([" ".join(col) for col in table_types])
    maf_create = "CREATE TABLE {}('index' INTEGER, {})".format(
        table_name, maf_table_def)
    con.execute(maf_create)
    # read and fill maf table
    maf_chunks = pd.read_table(maf_file,
                               dtype=maf_dtypes,
                               na_values=["."],
                               comment="#",
                               chunksize=chunksize)

    for chunk in maf_chunks:
        # process the chunk
        chunk['VAF'] = (chunk['t_alt_count'] /
                        (chunk['t_ref_count'] + chunk['t_alt_count']))
        if table_name in ('strelka', 'lancet', 'consensus'):
            chunk = chunk[common_cols]
        else:
            chunk = chunk[needed_cols]
        chunk.to_sql(table_name, con, if_exists='append')
    # create indexes
    print("Indexing table '{}'".format(table_name))
    for index_name, fields in maf_indexes:
        index_statement = "CREATE INDEX {table}_{name} ON {table}({field})".format(
            name=index_name,
            table=table_name,
            field=", ".join(fields)
        )
        con.execute(index_statement)

# Read the metadata file and load into a table called "samples"
if args.meta_file:
    print("Reading file {} to table 'samples'".format(args.meta_file))
    if args.overwrite:
        con.execute("DROP TABLE IF EXISTS samples")
    metadata_df = pd.read_table(args.meta_file)
    metadata_df.to_sql("samples", con)
