Known CNS signatures were obtained from the RefSig website: [https://signal.mutationalsignatures.com/explore/study/1
](https://signal.mutationalsignatures.com/explore/study/1
). The website no longer has a bulk download option, so the 8 mutational signatures were manually downloaded into separate csv files shown here, called `refsig_<NAME OF SIGNATURE>.csv`. The script `csv_to_rds_matrix.R` collects all CSV files into an appropriately-formatted matrix that can be used with `sigfit`.
