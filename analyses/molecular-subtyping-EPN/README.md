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

4. <b> `03-subgrouping_samples.py`  </b>  is a script that takes the table `results/EPN_all_data.tsv`  as input and adds a column that groups the samples into one of these groups - ST-EPN-RELA, ST-EPN-YAP1, PF-EPN-A, and PF-EPN-B. A new column named `subgroup` is added to the input table and saved in `resilts/EPN_all_data_withsubgroup.tsv`. It is possible that a sample is assigned more than one subgroup based on the conditions below. In those cases multiple subtypes are added to the subgroup column.  The logic for subtyping these are as follows - 
    - If fusion column values in [input](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/molecular-subtyping-EPN/results/EPN_all_data.tsv) table is higher  than 0, then the corresponding subtype is added to the last column. The  below table shows which fusion is associated with which Ependymoma subtype
            <table>
                <tr>
                    <th>Subtype name</th>
                    <th>Fusion genes</th>
                </tr>
                <tr>
                    <td>ST_EPN_RELA</td>
                    <td>C11orf95--RELA, LTBP3--RELA, PTEN--TAS2R1</td>
                </tr>
                <tr>
                    <td>ST_EPN_YAP1</td>
                    <td>YAP1--MAMLD1, C11orf95--MAML2, YAP1--FAM118B</td>
                </tr>
                <tr>
                    <td>PT_EPN_A</td>
                    <td>--</td>
                </tr>
                <tr>
                    <td>PT_EPN_B</td>
                    <td>--</td>
                </tr>
            </table>
     
    - If gene expression column for the genes given below under each subtype is greater than 3, then the subtype is added to the last column 
            <table>
                <tr>
                    <th>Subtype name</th>
                    <th>Gene expressions</th>
                </tr>
                <tr>
                    <td>ST_EPN_RELA</td>
                    <td>L1CAM, RELA</td>
                </tr>
                <tr>
                    <td>ST_EPN_YAP1</td>
                    <td>ARL4D, CLDN1</td>
                </tr>
                <tr>
                    <td>PT_EPN_A</td>
                    <td>CXorf67</td>
                </tr>
                <tr>
                    <td>PT_EPN_B</td>
                    <td>GPBP1</td>
                </tr>
            </table>
    -  If CNV columns have a value greater than 1, then the below subtype is assigned
            <table>
                <tr>
                    <th>Subtype name</th>
                    <th>CNV gain/loss</th>
                </tr>
                <tr>
                    <td>ST_EPN_RELA</td>
                    <td>9p_loss, 9q_loss</td>
                </tr>
                <tr>
                    <td>ST_EPN_YAP1</td>
                    <td>11q_loss, 11q_gain</td>
                </tr>
                <tr>
                    <td>PT_EPN_A</td>
                    <td>1q_loss</td>
                </tr>
                <tr>
                    <td>PT_EPN_B</td>
                    <td>6q_loss, 6p_loss</td>
                </tr>
            </table>  

        - ST_EPN_RELA subtype is associated with `CDKN2A loss`. If the column `consensus_focal_CN_CDKN2` shows loss or `GISTIC_focal_CN_CDKN2A` has a value less than 0.0, then those samples are associated with the ST_EPN_RELA subtype 

        -   The following formula was implemented for these two columns `breaks_density-chromosomal_instability_CNV` and `breaks_density-chromosomal_instability_SV` from input table and the values were added as column names `SV instability` and `CNV instability`

                `(break density value for sample - median) / interquartile range`   




