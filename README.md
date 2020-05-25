# Locational Choices
Models associated with the article: Blanchard, Simon J., Tatiana Dyachenko, Keri L. Kettle (2020). "Locational Choices: Modeling Consumer Preferences for Proximity to Others at Reserved-Seating Venues." <i>Journal of Marketing Research</i> [<A HREF="http://www.perceptionstudies.com/papers/Blanchard_PS_2020.pdf" target="_blank">PDF</A>].

## Code

In the folder <strong>PS</strong>, we provide the R codes of our Bayesian statistical model for the analysis of locational choice data. Two R functions and a DLL are provided. 

In the folder <strong>CNN_Benchmarks</strong>, we provide a naive implementation of covolutional neural networks (trained with adam) on the same data. The code was developed by <A HREF="https://www.gerad.ca/en/people/theo-moins">Theo Moins</A>. Taking a training file and holdout file as input, it outputs holdout predictive accuracy (top1, top5). 

## Datasets

The datasets were collected via <A HREF="http://www.seatmaplab.com" target="_blank">Seatmaplab.com</A>. 

<A HREF="http://www.seatmaplab.com/experiment/128" target="_blank">Try a study for yourself</A>.

For each set of analysis with CNNs, we generated a training (<b>in-sample</b>) and a <b>holdout</b> file. The characteristics of each dataset are described in more detail in Appendix A in the paper. Datasets, the raw outputs from seatmaplab and the R files to generate the training and holdout files can downloaded from the zip archives <A HREF="https://seatmaplab.com/public/locationalchoicedatasets/">here</A>.

## Usage

### PS (Bayesian Choice Models)

We provide three files: 
- FUN_PS_LPSreg_het_Cpp_withCov_diffNumCh.r
- FUN_PS_LPSreg_het_Cpp_withCov_diffNumCh_predictive.R
- PS_locational_wCov_noLambda_20180625.dll

<b>bayesm</b> and <b>Rcpp</b> are required. When using the functions, these three files are assumed to the in the working directory along with the <A HREF="https://seatmaplab.com/public/locationalchoicedatasets/">data</A> (see zip archives for sample usage and data processing).

### CNN_Benchmarks

We provide a Jupyter notebook. It requires Torch, numpy, csv, and pyplot. It requires a training and holdout csv file, available <A HREF="https://seatmaplab.com/public/locationalchoicedatasets/">here</A> (see CNN files column).
