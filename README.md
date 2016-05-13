# Foreign Fighter Supply

This repository has the replication files for the paper by George Derpanopoulos and Luke Sonnet on Foreign Fighter Supply and Kernel Regularized Hurdle Negative Binomial (KRHNB). It also contains the script to use KRHNB. Please try and replicate, and please tell us when it fails!

## Folder structure

* [`tex/`](tex/) contains the tex file, working paper, and bibliography
	* `tex/tabs_figs` contains tables, figures, and more produced by the R code for the .tex file
* [`code/`](code/) contains the analysis files, files that build the data, and some supporting functions
* [`krhnb/`](krhnb/) contains the script to run KRHNB
* [`savedata/`](savedata/) contains some saved .RData files because some of the analyses can take time
* [`data/`](data/) contains the cleaned data and the raw data used to create the cleaned data

## Main files

* [`code/krhnb.R`](code/krhnb.R) is the main script for the KRHNB method
* [`code/analyze_ff.R`](code/analyze_ff.R) is the main analysis script
* [`code/analyze_krhnb_performance.R`](code/analyze_krhnb_performance.R) is used to evaluate the OOS performance of our model
* [`data/foreignFightersImputed.csv`](data/foreignFightersImputed.csv) is the main imputed dataset
* [`tex/derpanopoulos_sonnet_ff.pdf`](tex/derpanopoulos_sonnet_ff.pdf) is the paper