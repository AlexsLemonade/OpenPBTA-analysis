# Molecular Subtyping of Ependymoma

<b>Module Authors:</b> Teja Koganti(<a href="https://github.com/tkoganti">@tkoganti</a>) and Josh Shapiro(<a href="https://github.com/jashapiro">@jashapiro</a>)

In this analysis we subtype ependymoma samples based on fusions, CNV, NFKB_pathway_GSEAscore, breaks_density-chromosomal_instability_CNV, breaks_density-chromosomal_instability_SV, GISTIC_focal_CN_CDKN2A and gene expression data

## Usage
`bash run-molecular-subtyping-EPN.sh`

This above  script is designed to change to this directory to run, so it should run from any location.

## Folder content

1. <b>`00-subset-for-EPN.R`</b> is a script that takes subsets only  ependymoma samples expression data for CI. The script uses `pbta-histologies.tsv` file to filter for ependymoma samples and     `pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds` file for expression  data

2. <b>`01-make_notebook_RNAandDNA.py`</b> is a script that filters WGS  and RNA-seq ependymoma samples from `pbta-histologies.tsv` file and adds `disease_group` column to the output file based on the primary_site. The values for `disease_group` are supratentorial/infratentorial/undetermined

3. <b>`02_ependymoma_generate_all_data.py`</b>  is a script that takes in expression, GISTIC, fusion, breakpoint, GISTIC, GSVA files to add values from these tables as new columns to the input notebook. Output from `01-make_notebook_RNAandDNA.py` script is used as input notebook. The output notebook from this is saved to `results/EPN_all_data.tsv`

4. <b> `03-subgrouping_samples.ipynb`  </b>  is a script that takes the table `results/EPN_all_data.tsv`  as input and adds a column that groups the samples into one of these groups - EPN, ST RELA, EPN, ST YAP1, EPN, PF A, and EPN, PF B. A new column named `subgroup` is added to the input table and saved in `results/EPN_all_data_withsubgroup.tsv`.
    - This script prioritizes features of subgroups first and does not assign those samples to any other subgroups. For example `RELA` fusions are prioritized for `EPN, ST RELA` subgroup and not assigedn to any other  groups. Two functions help achieve this in the script - 1) `prioritized_fusion` for `EPN, ST RELA` and `EPN, ST YAP1` groups and 2) `prioritizing_PT_EPN` for `EPN, PF A` and `EPN, PF B` groups

    - There is another function `subgroup_func` that checks that samples are not in `samples_assigned`  list  and checks the  value of certain columns and assigns subgroups to samples accordingly.

    - From the [input file here](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/molecular-subtyping-EPN/results/EPN_all_data.tsv) values for various columns are considered for assigning subgroups. Following are  the columns and values used (prioritized column shows `yes`, then samples  that have  that feature are _only_ assigned to that subgroup) -
            <table>
                <tr>
                    <th>Subtype name</th>
                    <th>Fusion genes</th>
                    <th>proiritized?</th>
                </tr>
                <tr>
                    <td>EPN, ST RELA</td>
                    <td>C11orf95--RELA, LTBP3--RELA, PTEN--TAS2R1</td>
                    <td>Yes</td>
                </tr>
                <tr>
                    <td>EPN, ST YAP1</td>
                    <td>YAP1--MAMLD1, C11orf95--MAML2, YAP1--FAM118B</td>
                    <td>Yes</td>
                </tr>
            </table>

    - Combination of gene expression and CNV gain/loss were checked below.
            <table>
                <tr>
                    <th>Subtype name</th>
                    <th>Gene expressions and CNV with  value</th>
                    <th>prioritized</th>
                </tr>
                <tr>
                    <td>EPN, PF A</td>
                    <td>CXorf67_expr_zscore>3</td>
                    <td>Yes</td>
                </tr>
                <tr>
                    <td>EPN, PF A</td>
                    <td>CXorf67_expr_zscore>3 and 1q_gain>0</td>
                    <td>Yes</td>
                </tr>
                <tr>
                    <td>EPN, PF A</td>
                    <td>TKTL1_expr_zscore>3 and 1q_gain>0</td>
                    <td>Yes</td>
                </tr>
                <tr>
                    <td>EPN, PF B</td>
                    <td>GPBP1_expr_zscore>3 and 6q_loss>0</td>
                    <td>Yes</td>
                </tr>
                <tr>
                    <td>EPN, PF B</td>
                    <td>GPBP1_expr_zscore>3 and 6p_loss>0</td>
                    <td>Yes</td>
                </tr>
                <tr>
                    <td>EPN, PF B</td>
                    <td>IFT46_expr_zscore>3 and 6q_loss>0</td>
                    <td>Yes</td>
                </tr>
                <tr>
                    <td>EPN, PF B</td>
                    <td>IFT46_expr_zscore>3 and 6p_loss>0</td>
                    <td>Yes</td>
                </tr>
            </table>
    -  Gene expressions and CNV values that  were considered for assigning subgroups. This  table contains featured that were _not_ prioritized meaning if a sample has `PTEN--TAS2R1` and `C11orf95--MAML2` fusions, that sample will show as subgrouped under both  `EPN, ST RELA`  and `EPN, ST YAP1`
            <table>
                <tr>
                    <th>Subtype name</th>
                    <th>CNV gain/loss</th>
                    <th>prioritized</th>
                </tr>
                <tr>
                    <td>EPN, ST RELA</td>
                    <td>PTEN--TAS2R1>0, <br/> 9p_loss>0, <br/> 9q_loss>0, <br/> RELA_expr_zscore>3, <br/> L1CAM_expr_zscore>3 </td>
                    <td>No</td>
                </tr>
                <tr>
                    <td>EPN, ST YAP1</td>
                    <td>C11orf95--MAML2>0, <br/> 11q_loss>0, <br/> 11q_gain>0, <br/> ARL4D_expr_zscore>3, <br/>CLDN1_expr_zscore>3</td>
                    <td>No</td>
                </tr>
            </table>  

      -   The following formula was implemented for these two columns `breaks_density-chromosomal_instability_CNV` and `breaks_density-chromosomal_instability_SV` from input table and the values were added as column names `SV instability` and `CNV instability`

                `(break density value for sample - median) / interquartile range`   
