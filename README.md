# MSB1015-Scientific-Programming-Gut-Microbiome-Research  
Period 1 encapsulating choosing a data set (gut microbiome) and conducting scientific programming analysis to answer research questions  

Research Question:  To what extent do longitudinal gut microbiome trajectories and inflammatory biomarkers predict individual dietary responses to a FV intervention, and what are the key microbial drivers that characterise the treatment group from the control group.

Data can be obtained through a publicly available website: https://datadryad.org/dataset/doi:10.5061/dryad.547d7wmfk#citations.  
Download the zip file -> unzip -> read to be used for your programming

Data contains: ASV excel file | Inflammatory Biomarkers csv file | Read me markdown file. 

## Programming Outline - The Steps Followed   

### Preprocessing  

1) data manipulation and organising - extracting relevant data and subset data into numerical/categorical values
2) data distribution - check and investigate the nature of the data, library depth, detected ASVs (Amplicon Sequence Variants)
3) prevalence and abundance distribution - discover prevalence and abundance distribution of data based on chosen thresholds
4) feature filtering - > 5% prevalence & > 0.01% abundance thresholds obtained for feature filtering and subset kept features
5) post filtering data representation - PCA to investigate data distribution and patterns
6) scree plots - showcases the PC (Principal Component) and PCo (Principal Coordinate) contributions for PCA (Principal Component Analysis) and PCoA (Principal Coordinate Analysis) respectively
7) sparsity profile and heteroscedasticity - validate and perform QC (Quality Control) to analyse the data behaviour - zero values and mean variance across abundance respectively.
8) rCLR normalisation - Apply robust Centred Log Ratio to the data to account for the volume of sparseness
9) PCoA - display the beta diversity and underlying feature patterns
10) Removing outliers - based on euclidean space centroids of distance metric with any point lying further than 3 SD (Standard Deviations) away

### Statistical Testing

11) PERMANOVA / R-PERMANOVA - testing differences in microbial composition between control and intervention groups <br />
R-PERMANOVA = Restricted Permutation Multivariate Analysis Of Variance
12) RM-ASCA - test for significance in microbial compositional shifts with respect to treatment and time <br />
RM-ASCA = Repeated Measure ANOVA Simultaneous Component Analysis

### Machine Learning

13) prepare the data - collect necessary data
14) (optional) parallel computing - more efficiently use your systems core to perform machine learning algorithms faster
15) imputation - random forest imputation method that learns on train-test split data to fill NAs
16) QC for imputation - density plots displaying that data is not distorted before/after imputation
17) apply group folds - Sample integrity needs to be preserved, treated as dependent blocks for each individual because each gave 6 samples which cannot be split (sample blocks become dependent on each other)
18) define data portioning - 80/20 data split for training and testing respectively
19) running machine learning models - EN, RF, and XGB (Elastic Net, Random Forest, XGBoost) respectively
20) model evaluations - extracting model accuracy, sensitivity, specificity, and AUC (Area Under Curve) of ROC (Receiver Operating Characteristic)
21) model results - compare each model and extract each model top contributing features

### Overview

22) extract key drivers - utilising RM-ASCA to get significant ASV species and to then check whether any predictive power within exists
23) check baseline (Week 0) and end-period (Week 8) microbiome composition differences between each sample iteration

### Conclusion & Limitations

24) XGB performed the best with regards to AUC
25) Specific significant features (ASVs) extract carry not much predictive power
26) A more concrete machine learning pipeline needs to be developed - since now this is more experimental
27) Deeper structure and further discoveries towards the nature of the data can impact machine learning performance
28) More developed model comparison metric is required
29) Experiment with not combining the ASV table and Inflammatory biomarker datasets for machine learning

## Software and Packages Required

### Software
Language: R (Version 4.6) <br />
Environment: RStudio <br />
https://docs.posit.co/ide/user/#rstudio-ide-oss-downloads?ref=rstudio-legacy

### Packages
All provided in the script and will be installed if the user does not have them, else it is skipped (all versions (V) additionally provided in the script) <br />
1) readxl - V 1.5.0.1
2) tidyverse - V 2.0.0
3) janitor - V 2.2.1
4) ggplot2 - V 4.0.3
5) BiocManager - V 3.23
6) pcaMethods - V 2.4.0
7) factoextra - V 2.2.0
8) compositions - V 2.0.9
9) vegan - V 2.7.6
10) caret - V 7.0.1
11) missForestPredict - V 1.0.1
12) doParallel - V 1.0.17
13) glmnet - V 5.1
14) ranger - V 0.18.0
15) xgboost - V 3.2.1.1
16) pROC - V 1.19.1
17) Ckmeans.1d.dp - V 4.3.6
18) randomForest - V 4.7.1.2
19) lme4 - V 2.0.6
20) lmerTest - V 3.2.1

## Running the script

The script is ready to be run without any errors or conflicts, assuming the necessary softwares and packages are installed <br />
The first line of code is where the directory of your unzipped files lie (the dataset) - therefore, place the file path in <br />
w_dir - working directory. <br />
Assuming the user has GitHub and uses git through the software RStudio, if you want the plots to be saved, then you must create a folder; <br />
use git command: mkdir -p plots <br />
Similarly it can be done to create other folders to save other objects, such as the source code <br />
<br /> 
If the user does not have git installed for their environment, then one can follow this link to set it up: https://jennybc.github.io/2014-05-12-ubc/ubc-r/session03_git.html <br />

Side Note: For any package that asks for updates when running the script, type "a" in the console, which updates "all" packages and you may proceed forward.




