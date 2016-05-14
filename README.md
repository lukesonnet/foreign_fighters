# Foreign Fighter Supply and Kernel Regularized Hurdle Negative Binomial

This repository has the replication files for the paper "Predicting Foreign Fighter Flows to Syria Using Machine Learning: An Introduction to Kernel Regularized Hurdle Negative Binomial" ([GitHub link](tex/derpanopoulos_sonnet_ff.pdf), [DB link](https://www.dropbox.com/s/1ny0cewhyv2o4vb/derpanopoulos_sonnet_ff.pdf?dl=0)) by George Derpanopoulos and [Luke Sonnet](http://lukesonnet.github.io/). It also contains the script to use KRHNB. Please try and replicate, and please tell us when it fails!

To fully replicate, you have to do the following:
* Make sure all of the packages we use are installed (may take some time to work through all of the files, we will look to automate and facilitate this in the future)
* Make sure you set the working directory appropriately in each of the `R` files
* Run [`code/build_ff_data.R`](code/build_ff_data.R) to build our full dataset from the raw data
* Run [`code/analyze_ff.R`](code/analyze_ff.R) to run the full analysis. In the script there will be some flags you will have to turn to `TRUE` in order for the full analysis to run; these processes or slow so we instead rely on cached data
* Run [`code/analyze_krhnb_performance.R`](code/analyze_ff.R) to run the OOS performance analysis. Again in the script there will be some flags you will have to turn to `TRUE` in order for the full tests to run; these processes or slow so we instead rely on cached data
* Run [`tex/derpanopoulos_sonnet.tex`](tex/derpanopoulos_sonnet.tex) to recreate our final PDF. Note that we edit the map to trim some white space and that we directly paste edited versions of the tables into the tex file. You will have to update those in the tex file to get them to update if you change some features of the code

# Folder structure

* [`tex/`](tex/) contains the tex file, working paper, and bibliography
	* [`tex/tabs_figs`](tex/tabs_figs) contains tables, figures, and more produced by the R code for the .tex file
* [`code/`](code/) contains the analysis files, files that build the data, and some supporting functions
* [`krhnb/`](krhnb/) contains the script to run KRHNB
* [`savedata/`](savedata/) contains some saved .RData files because some of the analyses can take time
* [`data/`](data/) contains the cleaned data and the raw data used to create the cleaned data

# Main files

* [`krhnb/krhnb.R`](code/krhnb.R) is the main script for the KRHNB method
* [`code/analyze_ff.R`](code/analyze_ff.R) is the main analysis script
* [`code/analyze_krhnb_performance.R`](code/analyze_krhnb_performance.R) is used to evaluate the OOS performance of our model
* [`data/foreignFightersImputed.csv`](data/foreignFightersImputed.csv) is the main imputed dataset
* [`tex/derpanopoulos_sonnet_ff.pdf`](tex/derpanopoulos_sonnet_ff.pdf) is the paper