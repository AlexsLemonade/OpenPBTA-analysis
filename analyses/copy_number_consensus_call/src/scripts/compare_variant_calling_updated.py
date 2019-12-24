################## PURPOSE #################

# This script is to do PAIR comparisons using the three 3 caller CNVs files
# This script takes in the 3 caller CNVs files and put their content in to a nested dictionary.
# There are 3 dictionaries within the bigger dictionary , each dictionary is the content of a singular caller method.

# For-loops are used to loop through these 3 files to do pair comparison of the 3 files.
# We compare manta - cnvkit
# then  manta - freec
# and finally cnvkit - freec
# We'll end up with 3 output files.

# The principle for finding consensus CNVs is that if any two CNVs overlap at least 50% reciprocally, we will
# take the common section of those overlapping CNVs as the new consensus CNV.

################# Output format ###############

#   column1    column2     column3     column4                     column5                     column6                  column7     column8     column9
#
#   chr#       consensus   consensus   raw manta coordinates       raw cnvkit coordinates      raw freec coordinates    cnv type    sample      filename
#              start_pos   end_pos     that made up this CNV       that made up this CNV       that made up this CNV                name

################# ASSUMPTION ###############

# None of the files have overlapping segments within its own file.

############################################




# Imports in the pep8 order https://www.python.org/dev/peps/pep-0008/#imports
# Standard library
import argparse
import sys
import os
import re

# Related third party
import numpy as np
import pandas as pd



########################## DEFINE FUNCTIONS ###############################

## Function to read content of input files
def read_input_file(input_file):
    ## Check to see if the CNVs bed file is empty, if so, make an empty variable
    if os.stat(input_file).st_size == 0:
        content_dict = {}

    ## If the file has content, then load in the content
    else:
        with open(input_file) as file:

            ## Strip the trailing spaces for each line
            stripped_content = [line.rstrip() for line in file.readlines()]

            ## Split each line up by any white space between columns
            fin_input_content = [re.split('\s+',line) for line in stripped_content]

            ## Create a list of unique keys/chromosome
            dict_key = []
            [dict_key.append(i[0]) for i in fin_input_content]
            dict_key_unique = np.unique(dict_key)

            ## Make the dictionary
            content_dict = {el:[] for el in dict_key_unique}

            ## Populate the dictionary
            {content_dict[line[0]].append(line[1:]) for line in fin_input_content}

    return content_dict

## Function to get consensus CNVs.
def generate_consensus(list1_name,list1,list2_name,list2):

    ## Define a variable to store the list of consensus CNVs
    consensus_list = []

    ## For every CNVs of list1, we look through all CNVs in list2 that overlaps and perform actions
    ## The use of dictionary is to speed this process up.\
    ## Loop through all chromosomes of list1
    for chr_list1 in list1:

        ## For every CNVs in that chromosome,
        for m in list1[chr_list1]:

            ## Assign the start_pos, end_pos, and copy number of list1 to variables
            start_list1 = int(m[0])
            end_list1 = int(m[1])
            cnv_list1 = m[2]

            ## Make sure end > start
            if start_list1 > end_list1:
                start_list1, end_list1 = end_list1, start_list1

            ## Assign variable for the CULMULATIVE coverage of the current CNV from list1
            coverage_list1 = 0

            ## Assign variables to store information for list2
            list2_overlap_len = 0
            list2_total_len = 0
            list2_start_list = []
            list2_end_list = []
            list2_chr_str_end = ''


            ## For every of the CNV in list one, we look over ALL CNVs in list 2 and make comparisons
            if chr_list1 in list2:
                for n in list2[chr_list1]:

                    ## Assign start_pos, end_pos, and copy number of list2 to variables
                    start_list2 = int(n[0])
                    end_list2 = int(n[1])
                    cnv_list2 = n[2]

                    ## Make sure end > start
                    if start_list2 > end_list2:
                        start_list2, end_list2 = end_list2, start_list2

                    ## For the current CNV in list1, if it overlaps with any CNV from list2
                    if start_list1 <= end_list2 and end_list1 >= start_list2:

                        ## Take the common region and store them in variables
                        start = max(start_list1,start_list2)
                        end = min(end_list1,end_list2)

                        ## For list1's CNV,
                        ## if there are any overlaps at all - immediately add the coverage into the overall coverage for that specific CNV of list1
                        coverage_list1 += (end - start + 1) / (end_list1 - start_list1 + 1)

                        ## For list2's CNV
                        ## If the overlapping covers >= 50% of its length,
                        ## then we add in the start, end coordinate, total overlap length, and total len to different lists
                        ## This is done to account for 1 CNV from list1 overlapping with MULTIPLE CNVs from list2
                        if (end - start +1) / (end_list2 - start_list2 + 1) >= 0.5:

                            ## Add the overlapping length to a list of overlap length
                            list2_overlap_len += (end - start + 1)

                            ## Add the total length to a list of total length
                            list2_total_len += (end_list2 - start_list2 + 1)

                            ## Add start_pos and end_pos to a list to choose from later on
                            list2_start_list += [start_list2]
                            list2_end_list += [end_list2]

                            ## Add the raw coordinates to a growing list that gets output later on to the final consensus file
                            list2_chr_str_end += str(cnv_list2)




                ## The coverage of list1 is already calculated as the "coverage_list1" variable
                ## ASSUMING there is no overlaping segments within 1 method
                ## then it is safe to add/accumulate the coverage_list1 value. Refer to assumptions (line 23)

                ## Next we calculate the coverage of list2 since there might be more than 1 overlapping CNV from list2
                if list2_total_len > 0:
                    coverage_list2 = list2_overlap_len / list2_total_len
                else:
                    coverage_list2 = 0

                ## Check to see if the overlap is >= 50% RECIPROCALLY
                if coverage_list1 >= 0.5 and coverage_list2 >= 0.5:
                    ## if it is 50% reciprocally, we chose from the list2 starts and ends
                    ## we take out the smallest start and biggest end
                    ## This maximize the CNV length from list2
                    fin_start_list2 = min(list2_start_list)
                    fin_end_list2 = max(list2_end_list)

                    ## Taking only the common region as the final segment.
                    ## Compare between list2 vs list1 start and end - take the common segment
                    chrom_start = max(fin_start_list2,start_list1)
                    chrom_end = min(fin_end_list2,end_list1)

                    ## Format the output consensus CNV
                    ## Column 4 is manta's raw coordinates
                    ## Column 5 is cnvkit's raw coordinates
                    ## Column 6 is freec's raw coordinates

                    ## if the files input are MANTA and CNVKIT, put info in column 4 and 5, 6th column is null
                    if (list1_name == 'manta' and list2_name == 'cnvkit') or (list1_name == 'cnvkit' and list2_name == 'manta'):
                        overlap_chrom = [chr_list1,str(chrom_start),str(chrom_end),str(cnv_list1).strip(','),list2_chr_str_end.strip(','),'NULL',m[-1]]

                    ## if the files input are MANTA and FREEC, put info in column 4 and 6, 5th column is null
                    elif (list1_name == 'manta' and list2_name == 'freec') or (list1_name == 'freec' and list2_name == 'manta'):
                        overlap_chrom = [chr_list1,str(chrom_start),str(chrom_end),str(cnv_list1).strip(','),'NULL',list2_chr_str_end.strip(','),m[-1]]

                    ## if the files input are CNVKIT and FREEC, put info in column 5 and 6, 4th column is null
                    elif (list1_name == 'cnvkit' and list2_name == 'freec') or (list1_name == 'freec' and list2_name == 'cnvkit'):
                        overlap_chrom = [chr_list1,str(chrom_start),str(chrom_end),'NULL',str(cnv_list1).strip(','),list2_chr_str_end.strip(','),m[-1]]

                    ## Add the formatted consensus CNV into a growing list of consensus CNVs
                    consensus_list = consensus_list + [overlap_chrom]

    return consensus_list


## Function to save output files
def save_to_file(output_file_content, output_path, sample_name):
    ## Open up the file to write in it, the 'w' option overwrites the file if the file exists.
    with open( output_path, 'w') as file:
        ## Loop through the output file and print line by line
        for k in output_file_content:

            ## Add the sample name and file name to each line
            single_name = output_path.split('/')[-1]
            file.write('\t'.join(k) + '\t'+ sample_name + '\t' + single_name +'\n')
    sys.stderr.write('$$$ Write to file ' + str(output_path) + ' was sucessful\n')



################### SCRIPT BODY #################################



## Define parser for the input file
parser = argparse.ArgumentParser(description="""This script takes in 3 bed files, each from one of the 3 callers
                                                 and find common CNVs between two files at a time.""")
parser.add_argument('--manta', required=True,
                    help='path to the manta file')
parser.add_argument('--cnvkit', required=True,
                    help='path to the cnvkit file')
parser.add_argument('--freec', required=True,
                    help='path to the freec file')
parser.add_argument('--manta_cnvkit', required=True,
                    help='path of the output consensus between manta and cnvkit')
parser.add_argument('--manta_freec', required=True,
                    help='path of the output consensus between manta and freec')
parser.add_argument('--cnvkit_freec', required=True,
                    help='path of the output consensus between cnvkit and freec')
parser.add_argument('--sample', default='no_sample',
                    help='sample name to use in the file')

args = parser.parse_args()


## Read in the input files as dictionaries
manta_dict = read_input_file(args.manta)
cnvkit_dict = read_input_file(args.cnvkit)
freec_dict = read_input_file(args.freec)

## Put the input dictionaries into a bigger dictionary
dict_of_input_content = {'manta':manta_dict, 'cnvkit':cnvkit_dict, 'freec':freec_dict}

input_file_names = list(dict_of_input_content.keys())

## Put the output file paths into their own lists that is to be iterated over
## This order is important
output_list = [args.manta_cnvkit, args.manta_freec, args.cnvkit_freec]

## Make a list of iteraby index for the MAIN for-loop below
## This give a list as followed: [0,1,2]
list_index = list(range(0,len(dict_of_input_content)))


## We have 3 items in dict_of_input_content. Each item is the content of an input file
## Make a list to store the final files
fin_list =[]

## Loop through list_index
for j, jval in enumerate(list_index):
    for k,kval in enumerate(list_index[j+1:]):

        ## Compare 1st and 2nd file from list_of_input_contents
        ## Then compare the 1st and 3rd
        ## Then compare the 2nd and 3rd
        list1 = dict_of_input_content[input_file_names[jval]]
        list2 = dict_of_input_content[input_file_names[kval]]

        ## Call the "generate_consensus" function to get consensus CNVs
        consensus_list = generate_consensus(input_file_names[jval],list1,input_file_names[kval],list2)



        ## Add the consensus list into a final list that hold the
        ## content of 3 files - manta_cnvkit , manta_freec, and cnvkit_freec
        fin_list.append(consensus_list)



## Print the output into files
for content, outfile in zip(fin_list, output_list):
    save_to_file(content, outfile, args.sample)

    ## Call the "save_to_file" function to print each consensus content to a file.

