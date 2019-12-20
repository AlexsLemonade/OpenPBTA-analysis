################## PURPOSE #################

# This script is to do PAIR comparisons using the three 3 caller CNVs files
# This script takes in the 3 caller CNVs files and put their content in to list.
# There are 3 elements in this list, each element is the content of a singular caller method.

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

## Put the input and output file paths into their own lists that is to be iterated over
## This order is important
list_of_files = [args.manta, args.cnvkit, args.freec]
list_of_output_files = [args.manta_cnvkit, args.manta_freec, args.cnvkit_freec]

## Create a list to store the content of the 3 caller CNV files
list_of_list = []


## Loop through the file paths and read in their contents to store in a list
for i in list_of_files:

    ## Check to see if the CNVs bed file is empty, if so, make an empty variable
    if os.stat(i).st_size == 0:
        fin_input_content = []

    ## If the file has content, then load in the content
    else:
        with open(i) as file:

            ## Strip the trailing spaces for each line
            stripped_content = [i.rstrip() for i in file.readlines()]

            ## Split each line up by any white space between columns
            fin_input_content = [re.split('\s+',i) for i in stripped_content]

            ## If the 1st column has '1' instead of 'chr1' for the chromosome number, add in 'chr'
            if fin_input_content[0][0].find('chr') == -1:
                for c,chromo in enumerate(fin_input_content):
                    fin_input_content[c][0] = 'chr' + fin_input_content[c][0]

    ## Add the content to a list
    list_of_list += [fin_input_content]


## Make a list of iteraby index for the MAIN for-loop below
## This give a list as followed: [0,1,2]
list_index = list(range(0,len(list_of_list)))


## We have 3 items in list_of_list. Each item is the content of an input file
## Make a list to store the final files
fin_list =[]

## Loop through list_index
for j, jval in enumerate(list_index):
    for k,kval in enumerate(list_index[j+1:]):

        ## Compare 1st and 2nd file from list_of_list
        ## Then compare the 1st and 3rd
        ## Then compare the 2nd and 3rd
        list1 = list_of_list[jval]
        list2 = list_of_list[kval]

        ## Define a variable to store the list of consensus CNVs
        consensus_list = []

        ## For every CNVs in list1, we look through the ENTIRETY of list2 and perform actions if 2 CNVs overlaps
        ## Could probably speed this up by using a dictionary.
        for m in list1:

            ## Assign the start_pos, end_pos, and copy number of list1 to variables
            chr_list1 = m[0]
            start_list1 = int(m[1])
            end_list1 = int(m[2])
            cnv_list1 = m[3]

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
            for n in list2:

                ## Assign start_pos, end_pos, and copy number of list2 to variables
                chr_list2 = n[0]
                start_list2 = int(n[1])
                end_list2 = int(n[2])
                cnv_list2 = n[3]

                ## Make sure end > start
                if start_list2 > end_list2:
                    start_list2, end_list2 = end_list2, start_list2

                ## For the current CNV in list1, if it overlaps with any CNV from list2
                if m[0] == n[0] and start_list1 <= end_list2 and end_list1 >= start_list2:

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

                ## if jval == 0 (manta) and kval == 1 (cnvkit), put info in column 4 and 5, 6th column is null
                if jval == 0 and kval == 1:
                    overlap_chrom = [m[0],str(chrom_start),str(chrom_end),str(cnv_list1).strip(','),list2_chr_str_end.strip(','),'NULL',m[-1]]

                ## if jval == 0 (manta) and kval == 2 (freec), put info in column 4 and 6, 5th column is null
                elif jval == 0 and kval == 2:
                    overlap_chrom = [m[0],str(chrom_start),str(chrom_end),str(cnv_list1).strip(','),'NULL',list2_chr_str_end.strip(','),m[-1]]

                ## if jval == 1 (cnvkit) and kval == 2 (freec), put info in column 5 and 6, 4th column is null
                elif jval == 1 and kval == 2:
                    overlap_chrom = [m[0],str(chrom_start),str(chrom_end),'NULL',str(cnv_list1).strip(','),list2_chr_str_end.strip(','),m[-1]]

                ## Add the formatted consensus CNV into a growing list of consensus CNVs
                consensus_list = consensus_list + [overlap_chrom]

        ## Add the consensus list into a final list that hold the
        ## content of 3 files - manta_cnvkit , manta_freec, and cnvkit_freec
        fin_list += [consensus_list]

## Print the output into files
for i,ival in enumerate(fin_list):
    try:
        ## Attempt to create a the file, if the file already exists, throw ERROR
        with open( list_of_output_files[i] , 'x') as file:
            sys.stderr.write('$$$ Created new file succesfully\n')

        ## Open up the file to write in it
        with open( list_of_output_files[i], 'w') as file:
            ## Loop through the output file and print line by line
            for k in ival:

                ## Add the sample name and file name to each line
                single_name = list_of_output_files[i].split('/')[-1]
                sample_name = args.sample
                file.write('\t'.join(k) + '\t'+ sample_name + '\t' + single_name +'\n')
        sys.stderr.write('$$$ Write to file ' + str(list_of_output_files[i]) + ' was sucessful\n')
    except:
        sys.stderr.write('!!! file ' + str(list_of_output_files[i]) + ' seems to be already existed OR the output path doesn\'t exist. Please check and try again\n')
        sys.exit()
