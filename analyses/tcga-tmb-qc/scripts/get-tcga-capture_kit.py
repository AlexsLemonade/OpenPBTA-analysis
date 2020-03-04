################## PURPOSE #################
#
# retrieve the TCGA's exome capture kit from GDC's file API endpoint 
# so that we can use the correct BED for process like somatic calling 
# and TMB calculation or other related analysis
#
############################################

import requests
import json
import pandas as pd
import os

# 1. get TCGA manifest from the data release
## set data bucket base url and version 
PBTA_BUCKET = 'https://s3.amazonaws.com/kf-openaccess-us-east-1-prd-pbta/data/'
RELEASE = 'release-v14-20200203/'

## GET and load manifest
tcga_manifest = requests.get(
    PBTA_BUCKET+RELEASE+"pbta-tcga-manifest.tsv").content

## iterate TCGA manifest, to get all the file names
tcga_manifest_lines = tcga_manifest.split("\n")
tcga_filenames = []
for line in tcga_manifest_lines:
    tcga_filenames.append(line.split("\t")[0])

# 2. hit GDC file API endpoint to get the details of the capture kit
## set GDC API base url and request headers
gdc_url = 'https://api.gdc.cancer.gov/files'
headers = {'Content-Type': 'application/json'}

## API request field, removed "analysis.metadata.read_groups.read_group_name"
## can add that back for details
fields = [
    'file_name',
    'analysis.metadata.read_groups.target_capture_kit_name',
    'analysis.metadata.read_groups.target_capture_kit_target_region'
]
fields = ','.join(fields)

## API request body 
payload = {
        'filters':{
            'op':'=',
            'content':{
                'field':'file_name',
                'value':tcga_filenames}},
        'format':'json',
        'fields':fields,
        'size':5000 # make sure we get all the returns
}
payload = json.dumps(payload)

## hit GDC API file endpoint
gdc_response = requests.post(gdc_url, headers=headers, data=payload)

# 3. handle GDC API return to find out capture kit url
gdc_response = gdc_response.json()
capture_kits = []

## iterate .data.hits entity manifest
for i in gdc_response['data']['hits']:
    for j in i['analysis']['metadata']['read_groups']:
        capture_kits.append([
            i['file_name'], 
            j['target_capture_kit_name'],
            j['target_capture_kit_target_region']
        ])

# 4. load capture kit into data frame, find unique kit download url
df = pd.DataFrame(capture_kits).drop_duplicates()
df.columns = ['filename','kit_name','kit_url']

# 5. output the capture kit data frame
df.to_csv(
    os.path.join('results', 'tcga-capture_kit-info.tsv'),
    sep='\t',
    index=False
)
