

wget https://gdc.cancer.gov/system/files/authenticated%20user/0/dtt-ui_v0.5.4_Ubuntu_x64.zip

unzip gdc-client_v1.0.1_Ubuntu14.04_x64.zip ./gdc-client

cp -pi ./gdc-client /usr/local/bin (if this does not work)
sudo cp -pi ./gdc-client /usr/local/bin

gdc-client download -m tcga_data/gdc_manifest_maf_files.txt

