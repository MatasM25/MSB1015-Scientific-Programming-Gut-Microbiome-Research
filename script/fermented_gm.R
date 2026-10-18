# INFO - R version: 4.6.0

# Link to the original research paper and data information
# DOI: https://doi.org/10.5061/dryad.547d7wmfk

# Set working directory for easy navigation
w_dir <- setwd("/Users/matasm/Downloads/MSc System Bio/Maastricht Uni/Y2/P1/Scientific Programming/doi_10_5061_dryad_547d7wmfk__v20241009")
dir() # proves you have all files necessary

# load data and explore
# README file and the website have full description of dataset features

install.packages("readxl", type = "source") # install and load packages to read excel file;
library(readxl)
packageVersion("readxl") # the version of package used for this project: 1.5.0.1

#
# NOTE: for all packages, if requested you may update all/some/none, does not matter what is chosen, in this project "a" is selected
#
# PLOTS are labelled in a trivial numerical order just for simplicity in saving them to a folder

seq_data <- read_xlsx("ASV_table_dryad.xlsx") 

# output of NGS via Illumina MiSeq -> microbial alpha, 
# beta diversity, and amplicon sequencing variants 
# (ASVs = exact sequence variants, inferred single DNA sequences from high-throughput marker gene analyses)

inf_markers <- read.csv("Inflammatory_markers.csv", na = c(".", "", "NA"))
# DRYAD specifies that "." represents an NA value = missing

# inflammatory markers measured in serum collected at two time points: at baseline and at 8 weeks
# 1 = week 0 | 2 = 8 weeks [in the datasets]


# RQ - Can week 0 (baseline) gut microbiome composition of patients predict heterogeneous individual responses to fermented vegetable diets at 8 weeks.
# Evolution:
# RQ - To what extent do longitudinal gut microbiome trajectories and inflammatory biomarkers predict individual dietary responses to a FV intervention, and what are the key microbial drivers that characterise the treatment group from the control group.


# 1 - PREPROCESSING DATA
#========================

# Since the data imported has no description from the website/article if the data is normalised, it needs to be checked.
# Normalising data is required because we cannot directly see which species directly affects the results, due to each species of the whole taxonomic composition is not equal.
# Therefore, the effect of each species for later research is not "directly" measurable/compared.

# NOTE: careful with analysis as there are different time points in the dataset;
# So one cannot use PERMANOVA or a random test-train split to feed into ML algorithms => data leakage

# raw data vs transformed;

install.packages("tidyverse") # data manipulation package
library(tidyverse)
packageVersion("tidyverse") # Version: 2.0.0

if(!require("janitor")) install.packages("janitor") # data cleaning package
library(janitor)
packageVersion("janitor") # Version: 2.2.1

# Relabel the variables for easier readability;
seq_data <- seq_data %>%
  mutate(time = factor(time, levels = c(1, 2), labels = c("W0", "W8")), Treatment = factor(Treatment, levels = c(1, 2), labels = c("UD", "FV"))) 
# relabel 1 = Week 0 & 2 = Week 8
# relabel 1 = Usual Diet & 2 = Fermented Vegetables

# The article website mentions that the first 18 columns is everything other than ASVs (Amplicon Sequencing Variants)
# Split to meta and ASVs for easier interpretability (not necessary)

seq_meta <- seq_data[1:18]
seq_asv <- seq_data[19:ncol(seq_data)]
dim(seq_meta) + dim(seq_asv) == dim(seq_data) # TRUE shows the columns have the same number = no data lost
# FALSE is rows + rows is 2 x rows

# When checking how many unqiue IDs there are there should be 87 as mentioned in the paper, but only 81 is given?
# n_distinct(seq_data$ID)

# Check how many samples of stool are for each individual;

samples_per_indv <- seq_data %>%
  count(ID, time, name = "n_samples")

# visualise data to see the distribution of samples to later utilise it for correct species contribution from individuals.
# similarly to variation on microbiome compositions fractions

# QUALITY CONTROL
#================

if(!require("ggplot2")) install.packages("ggplot2") # plotting package
library(ggplot2)
packageVersion("ggplot2") # Version: 4.0.3

p1 <- ggplot(samples_per_indv, aes(x = n_samples, fill = time)) + geom_bar() + facet_wrap(~ time) + labs(title = "Number of microbiome samples per patient", x = "Samples", y = "Number of patients")
ggsave("plots/p1.png")

# displays completely equal distribution from baseline to the end of the study.

# Now inspecting raw ASVs and to make sure they are numeric values, else any matrix operations will fail;

seq_asv_mat <- seq_asv %>%
  mutate(across(everything(), as.numeric))

# There should be no negative values, since this is microbiome species abundance

# Additionally, it is good to visualise library depth, it shall show which species contribute the most and leaste, in case some library biases need to be accounted for

seq_lib_sizes <- rowSums(seq_asv_mat, na.rm = TRUE)

p2 <- ggplot(data.frame(seq_lib_size = seq_lib_sizes), aes(x = seq_lib_size)) + geom_histogram(bins = 35) + labs(title = "Raw sequencing depth", x = "Total ASV count per sample", y = "Number of samples")
ggsave("plots/p2.png")
# The sequencing depth follows a close path to a normal distribution meaning, usually poisson or negative binomial distribution is expected
# Due to microbiome species represent their own discrete, independent activities

# This could mean that the data is already normalised. 
# Another check list is that the gut microbiome compositions are in percentage contributions that sunm to 1 that can be seen in "seq_lib_sizes"

# Locate the species that have no contribution to the total microbiotic composition.
# Depending on the fraction of "inactive" species, will act as an indicator for filtering;

# looking at the number of detected ASVs per sample to see what the contributes the most in the samples
detected_asv <- rowSums(seq_asv > 0)
summary(detected_asv) # shows the summary of ASVs in samples

p3 <- ggplot(data.frame(detected_asv = detected_asv), aes(x = detected_asv)) + geom_histogram(bins = 35) + labs(title = "Number of detected ASVs per sample", x = "Number of ASVs detected", y = "Number of samples")
ggsave("plots/p3.png")
# the distribution is not so skewed to any given direction. Follows an imperfect normal distribution, where 182 ASVs detected as shown in the summary is the highest for the frequency of samples

#^^^^^^^^^^^^^
# Exploratory
#^^^^^^^^^^^^^


#===================


#===============
# PREPROCESSING
#===============


# From the ASV tables, we have split the meta and ASV features from the table in two:
# seq_meta and seq_asv respectively

# We need to check the data for sparse values, sample distribution and abundances, to understand what samples may be hindering biological signal

#========================
# UNDERSTANDING THE DATA
#========================

# This is just as discussed in Week 2 journal club that we need to perform the initial data diagnosis;
# Check the direct problematic issue if any;

dim(seq_asv_mat) # total features from original table = 3063, where 18 features are meta => 3063 - 18 = 3045
# seq_asv_mat = seq_asv as only a numeric function was applied and luckily all data for asv was numneric
sum(is.na(seq_asv_mat)) # 0 means there are no missing values - NAs
sum(as.matrix(seq_asv_mat) < 0) # 0 means there are no negative values - so no negatively sequenced ASVs
sum(as.matrix(seq_asv_mat) > 1) # Similarly no value greater than 1 -> data is normalised in some ratio form
row_sums <- rowSums(seq_asv_mat) # returns 1 => definitely normalised and all features in a sample contribute some percentage

p4 <- ggplot(data.frame(row_sum = row_sums), aes(x = row_sum)) + geom_histogram(bins = 30) + labs(title = "Distribution of ASV relative abundance", x = "Sum of ASV abundances per sample", y = "Number of samples")
ggsave("plots/p4.png")
#distribution is relatively normal, but we see some blips which can be investigated, but also extremely low abundance that can be filtered out downstream

# I want to create a summary table of statics for my data;
# Looking into distribution, abundance, depth, and sparsity;

sample_qc <- seq_meta %>%
  mutate(seq_depth = rowSums(seq_asv_mat), detected_asvs = rowSums(seq_asv_mat > 0), zero_frac = rowMeans(seq_asv_mat == 0))
# The entries for detected_asvs show the distribution of sample features (ASVs) from the Illumina sequencing device

# These summaries just back up the claims mentioned before
summary(sample_qc$seq_depth)
summary(sample_qc$detected_asvs)
summary(sample_qc$zero_frac) # This identifies what the distribution of undetected samples

# I do not see the point in looking further into the sequencing depth since the data is showing relative abundance, so it is not raw data.
# Meaning, the plots are not showing what we want directly

p5 <- ggplot(sample_qc, aes(x = detected_asvs)) + geom_histogram(bins = 30) + labs(title = "Number of detected ASVs/sample", x = "Detected ASVs", y = "Samples")
ggsave("plots/p5.png")
# Displays 1 main peak, which can be seen from the summary at 182 detected ASVs following a normal distribution, and another smaller peak around 250 detected ASVs
# Additionally, showing some "low-quality" samples that need to be filtered with little ASBV abundance in samples

# Interestingly checking if the detected ASVs differ in the control and treated group;

p6 <- ggplot(sample_qc, aes(x = Treatment, y = detected_asvs, fill = Treatment)) + geom_boxplot() + labs(title = "Distribution of detected ASVs for UD vs FV", y = "Number of detected ASVs")
ggsave("plots/p6.png")
# While it is noticeable that the median distribution of ASVs in both groups is very similar, the direct spread is from the lower quartile range is different.
# The total spread from the whiskers is also larger for the control group, biologically this could mean;
# Since there is no dietary restrictions, the gut microbiome has a higher individual variability in microbial richness - seeing both very low and very high microbial richness.
# Whereas, the FV group could be exposed to much more specific microbiome species with fermented vegetables

# After some initial conversations, microbiome data needs to be checked for sparsity.
# Some ideas are to check the general spread of 0 abundant species in samples and create a threshold i.e. if there are more than 60% 0 abundant taxa than the feature must be dropped;

# So first I want to check the general prevalence of ASVs in samples [complete the average of each feature];
seq_asv_prev <- colMeans(seq_asv_mat > 0)
# and combine the data via making a data frame with the names of species and there mean prevalence;
prev_df <- tibble(ASV = names(seq_asv_prev), prevalence = seq_asv_prev)
# Opening the table, instantly we see ASV1 and ASV2 have ~90% and 100% prevalence in the data for samples, so they are always around

# We can conduct analysis for the actual total average of ASV abundance by including 0 values;
prev_df <- prev_df %>%
  mutate(avg_abundance = colMeans(seq_asv_mat), max_abundance = apply(seq_asv_mat, 2, max), zero_abundance = nrow(seq_asv_mat) - prevalence * nrow(seq_asv_mat))
# Abundance vs Prevalence;
# AB - measures total number of species present
# PR - measures the proportion of species present
# Now we have both mutated into the data frame
# We also look for the maximal abundance via "apply" where 2 stands for column (row is 1), so it is sifting the maxima in all columns and returning into the dataframe
# Lastly, looking for the amount of zero values in our matrix;
# by taking the row number and subtracting the row number multiplied with the prevalence proportion gives the result of species that are not "prevalent" - zero count

# Now we have a data frame where we can compare the entries;
# AB vs PR;
p7 <- ggplot(prev_df, aes(x = prevalence, y = avg_abundance)) + geom_point(alpha = 0.3, colour = "red") + scale_y_log10() + labs(title = "ASV prevalence vs abundace before filtering", x = "Prevalence", y = "Mean abundance - relative")
ggsave("plots/p7.png")
# I added the log scale to see the spread better, else it is very clumped together...
# What this displays: For low values of prevalence, we have "rare" species of ASVs and much more rare species than not, since the concentration lies in the lower side of prevalence
# Essentially we are seeing a logarithmic behaviour for rarity of species decreasing when increasing on the abundance scale

# We can also look at the prevalence distribution of ASVs;
p8 <- ggplot(prev_df, aes(x = prevalence)) + geom_histogram(bins = 30) + labs(title = "Distribution of ASV prevalence (rarity)", x = "Proportion of samples containing ASVs", y = "ASV count")
ggsave("plots/p8.png")
p9 <- ggplot(prev_df, aes(x = prevalence)) + geom_histogram(bins = 30) + scale_y_log10() + labs(title = "Distribution of ASV prevalence (rarity)", x = "Proportion of samples containing ASVs", y = "ASV count (log scale)")
ggsave("plots/p9.png")
# Original vs log scaled comparison, I re-plotted as "zoomed-in" version to see the distribution better
# The fraction of samples containing less ASVs is lower in numbers compared to the fraction of samples containing more ASVs

# I can now deduct that the data tells me that there more than 1000s of ASVs are present in ~5% of samples, showing the division in core microbiome composition and the transient microbiome
# On the right side we have the extremely prevalent species that are labelled "core" microbiomes for example microbiome species in charge of breaking down Short Chain Fatty Acids (SCFAs) that every human needs

# But I need to be careful, because batch effect or other confounding instances have not been accounted, there could also be signals showing noise, so we still need to complete these checks

# To move onto filtering, it is not as simple as expected;
# Since just choosing a threshold to just remove samples is usually the easiest, but for example 10% prevalence threshold would remove a 9% prevalent feature, which could be very much an important biological signal for FV or even UD groups
# Given some information at hand, I will use a 5% threshold for prevalence and 0.01% for abundance in filtering;

# We first require a data frame of a numeric matrix form;
seq_asv_df <- as.data.frame(seq_asv_mat)
# We want to extract the mean prevalence of species across samples, so we take the means from each column in the data frame;
prev_all <- colMeans(seq_asv_df > 0) # (GLOBAL) for all features that are non-zero, we will leave this for sparsity
prev_grouped <- seq_asv_df %>%
  mutate(Treatment = seq_meta$Treatment) %>% # Add the treatment column from meta data table to newly created prev_grouped
  group_by(Treatment) %>% # group by either FV or UD
  summarise(across(where(is.numeric), ~ mean(.x > 0))) %>% # return the prevalence column mean statistics for each ASV and treatment or control group
  select(-Treatment)
# is where(is.numeric) is necessary here because it will try take the mean of categorical results in the Treatment column

# Opening the new data frame we see for example the first entry of ASV1 is present in ~97.5% control group (UD) samples
# Secondly, ASV1 is present in ~ 86% of treatment group (FV) samples
# Explaining the mean prevalence

# Now we can selected which samples we want to keep;
keep_prev <- apply(prev_grouped >= 0.05, 2, any) # selecting value within this threshold from columns (2) in the data frame
keep_prev_test <- (prev_all >= 0.05) | apply(prev_grouped >= 0.05, 2, any) # with this OR operator - keep the prevalence of all global ASVs or ASVs grouped by treatment group within the threshold

table(keep_prev == keep_prev_test) # just a sanity check that both the outputs are the same: TRUE = 3045 

# for the moment the output dimensions for keep_prev is the same as prev_grouped, which is not what I want because I am trying to hand select ASVs that pass the threshold.
# A potential fix could be that biological noise is inflating the evaluation and we need to insert a detection limit i.e. not > 0, but > 0.00001 or so

# This was a false alarm, of course the dimension is the same, one can check if the values have been assessed by entering in the console: table(keep_prev) to check the count of TRUE/FALSE
# The only point that raises my eyebrow now, is that based on this threshold I would have a lot of samples removed;

sum(keep_prev == FALSE) / length(keep_prev) * 100 # shows more just over 73% of ASVs would be removed due to the 5% threshold, is this statistically and biologically viable?

# Similarly for abundance values;
abundance_all <- colMeans(seq_asv_df)

abundance_group <- seq_asv_df %>%
  mutate(Treatment = seq_meta$Treatment) %>%
  group_by(Treatment) %>%
  summarise(across(where(is.numeric), mean)) %>% # we are only calculating the general mean for abundance of sample species
  select(-Treatment)

keep_abundance <- apply(abundance_group >= 0.0001, 2, any)

sum(keep_abundance == FALSE) / length(keep_abundance) * 100 # Moreover, for the abundance a similar story, 79% of samples would have to be removed...

# Based on this post: https://www.biostars.org/p/9615831/ it is a reasonable choice of thresholds for prevalence and abundance.
# This paper (https://pmc.ncbi.nlm.nih.gov/articles/PMC10291178/#s2) on low abundance OTU filtering methods on reliability and variability of microbial composition experienced removal of around 76% OTUs, which is similar to my results here.

# Now we want to combine from the filtered data what is left; (prevalent ASVs are also abundant ASVs)
keep_asv <- keep_prev & keep_abundance
# Output the filtered data;
seq_asv_filtered <- seq_asv_mat[, keep_asv] # what this does is takes the matrix seq_asv_mat values and subsets them together with the columns of keep_asv, hence being in the second position of the square bracket

# Now we have our filtered ASVs based on prevalence and abundance thresholds
# Just for clarity I want to duplicate this filtered matrix (not necessary), only that it now contains the sample IDs;
seq_asv_filtered_ID <- seq_asv_filtered %>% # I can just add the sample ID column from meta data to this filtered matrix, since no samples where lost;
  mutate(ID = seq_meta$ID) %>% # add the ID from meta data to filtered data
  select(ID, everything()) # select the ID column to be represented at the start, not last - it is not possible to have them as row names since the sample IDs are non-unique, they are replicates of each other since multiple stool samples where collected from the same individual

# Repeating previous step for the filtered data prevalence and abundance distribution;
seq_asv_f_prev <- colMeans(seq_asv_filtered > 0)
f_prev_df <- tibble(ASV = names(seq_asv_f_prev), prevalence = seq_asv_f_prev)
f_prev_df <- f_prev_df %>%
  mutate(avg_abundance = colMeans(seq_asv_filtered))

p10 <- ggplot(f_prev_df, aes(x = prevalence, y = avg_abundance)) + geom_point(alpha = 0.3, colour = "red") + scale_y_log10() + labs(title = "ASV prevalence vs abundace after filtering", x = "Prevalence", y = "Mean abundance - relative")
ggsave("plots/p10.png")
# This displays the difference the filtering has on the data;
# We have more of a spread of features across prevalence and abundance scales. We have removed the dense 0-valued cluster and we have a more statistically relevant spread of rare and non-rare species.

# Now we can visualise the distribution of our filtered data to be more familiar with our results;
# For microbiome data, PCA is a good start, but PCoA is mostly used and I have used it before that can be used after transformation;

# By research and previous knowledge, CLR transformation - Centered Log Ratio will help remove sparsity;

# Having the filtered data from prevalence and abundance from microbiome features and taking the union of the two filters as the single dataframe, and make it numeric to work with a matrix;
seq_asv_f_mat <- as.matrix(seq_asv_filtered)

# checks if your system has BiocManager for the use of Bioconductor package installations, some packages throw errors even with the correct biocManager version, hence here just for consistency
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
update.packages(ask = FALSE, type = "source")
BiocManager::install(version = "3.23", type = "source")


# Packages necessary for PCA plotting;
BiocManager::install("pcaMethods", force = TRUE)
library(pcaMethods)
packageVersion("pcaMethods") # Version: 2.4.0

seq_asv_f_pc <- pca(seq_asv_f_mat, nPcs = 10) # arbitrarily choosing 10 principle components, since I am just checking the data behaviour
# This calculates the PCs for plotting, variance, scores, R2;
treatment_colour <- ifelse(seq_meta$Treatment == "UD", "red", "blue") # giving colour to the 2 treatment groups for nicer PCA visual;
plotPcs(seq_asv_f_pc, c(1, 2), col = treatment_colour) # this provides me with information of PC variance and relationship/distribution of samples within the chosen PC 1 and 2 in this case

# We can plot the scores also to see a more direct distribution;
asv_pca_scores <- as.data.frame(scores(seq_asv_f_pc))
p11 <- ggplot(asv_pca_scores, aes(x = PC1, y = PC2, colour = seq_meta$Treatment)) + geom_point(size = 1) + labs(title = "PCA plot of PC1 vs PC2 coloured by treatment") + stat_ellipse(level = 0.95) + scale_color_discrete(name = "Treatment")
ggsave("plots/p11.png")
# This gives us the data distribution and we can see the huge spread and no common pattern emerging.
# Additionally, from a visual basis there seems to be no global anomaly. Data points are heavily clustered around the origin, but not perfectly.
# I can add a ellipses to statistically see the points density behaviour, with 95% of the data being encapsulated.

#---------------------------
# Principal Component check
#---------------------------

# check where the variance lies in this analysis, to see which Principal Components (PCs) contribute the most - scree plot
if(!require("factoextra")) install.packages("factoextra") # plotting package via eigenvalues
library(factoextra)
packageVersion("factoextra") # Version: 2.2.0

# factoextra expects a prcomp - class of a data object, so we have to extract the standard deviation and set it as a prcomp variables, else fviz_eig cannot process the pca() object
mock_pc <- list(sdev = sDev(seq_asv_f_pc))
class(mock_pc) <- "prcomp"
# factoextra package function that allows to use the eigenvalues to plot the percentage of variance explained with labels for the proposed data;
p12 <- fviz_eig(mock_pc, addlabels = TRUE, barfill = "blue", barcolor = "blue", linecolor = "red", main = "Scree plot showcasing the PC's variance explained", ylim = c(0, 30)) # added a y limit because by default the percentages and the title overlap :(
ggsave("plots/p12.png")
# PC1 vs 2 with variance explained labels;
pc_variance <- (mock_pc$sdev^2 / sum(mock_pc$sdev^2)) * 100 # collecting the cummulative sum for variance explained across the PCs
p13 <- ggplot(asv_pca_scores, aes(x = PC1, y = PC2, colour = seq_meta$Treatment)) + geom_point(size = 1) + labs(title = "PCA plot of PC1 vs PC2 coloured by treatment", x = paste0("PC1 (", round(pc_variance[1], 1), "%)"), y = paste0("PC2 (", round(pc_variance[2], 1), "%)")) + stat_ellipse(level = 0.95) + scale_color_discrete(name = "Treatment")
ggsave("plots/p13.png")
# We see that the biggest contributors are PC1, 2, and 3. Therefore, we can also plot the relationship between PC1 and 3;
p14 <- ggplot(asv_pca_scores, aes(x = PC1, y = PC3, colour = seq_meta$Treatment)) + geom_point(size = 1) + labs(title = "PCA plot of PC1 vs PC3 coloured by treatment", x = paste0("PC1 (", round(pc_variance[1], 1), "%)"), y = paste0("PC3 (", round(pc_variance[3], 1), "%)")) + stat_ellipse(level = 0.95) + scale_color_discrete(name = "Treatment")
ggsave("plots/p14.png")
# relationship between PC2 and 3; (variance explained by PC2 and 3 are very similar)
p15 <- ggplot(asv_pca_scores, aes(x = PC2, y = PC3, colour = seq_meta$Treatment)) + geom_point(size = 1) + labs(title = "PCA plot of PC2 vs PC3 coloured by treatment", x = paste0("PC2 (", round(pc_variance[2], 1), "%)"), y = paste0("PC3 (", round(pc_variance[3], 1), "%)")) + stat_ellipse(level = 0.95) + scale_color_discrete(name = "Treatment")
ggsave("plots/p15.png")
# Most noticeably, the data is not very mean centered and the data is heavily distributed in "random" directions, which can be negated with CLR transformation

# Before any transformations, we need to account for the sparsity and fix it accordingly;
qc_metrics <- data.frame(Mean = colMeans(seq_asv_f_mat), Variance = apply(seq_asv_f_mat, 2, var), Zero_frac = colMeans(seq_asv_f_mat == 0)) %>% # all statistics for sparsity selection, feature names, column means, variance of features, averages of 0 values in matrix
  mutate(CV = sqrt(Variance) / Mean, Log_Mean = log10(Mean + 1), Log_Mean_p = log10(Mean + 10^-5)) %>% # Coefficient of variation, log mean is scaling the mean to account for 0 values with a pseudo-count to not compute log(0) - non-existent
  arrange(desc(Mean)) %>% # arrange the mean values from largest to smallest - descending order
  mutate(Cumulative_Reads = cumsum(Mean) / sum(Mean)) # calculates the running sum of the Mean vectors, this can help highlight the sequencing depth efficiency and assess microbiome diversity profiles

# Next, we can plot abundance to see the available mean ASV distribution against microbiome zero-reads;
p16 <- ggplot(qc_metrics, aes(x = Log_Mean, y = Zero_frac)) + geom_point(alpha = 0.4, colour = "red") + geom_smooth(method = "loess", colour = "green", se = FALSE) + labs(title = "Feature species Sparsity Profile", x = "Log10(Mean abundance + 1)", y = "Fraction of zero values - sparsity")
ggsave("plots/p16.png")
# https://r-statistics.co/geom_smooth-in-R.html - this page explains "loess" - automatic flexible local polynomial finder -> more for exploratory purposes

# Displayed on the plot; a lot of sparsity is located in the top left corner, so for low mean abundance we have high density of sparseness still occurring even after the prevalence and abundance filtering
# Highly abundant feature (to the right) have low zero fraction, so are picked up by the reads consistently. While the opposite is seen for very low mean abundance of features (top left).
# This could mean that many reads are zero-infalted or where not picked up by the Illumina sequencing device

# This is a more accurate plot of the previous version on a more biologically and mathematically relevant scale, showing true data behaviour;
p17 <- ggplot(qc_metrics, aes(x = Log_Mean_p, y = Zero_frac)) + geom_point(alpha = 0.4, colour = "red") + geom_smooth(method = "loess", colour = "green", se = FALSE) + labs(title = "Feature species Sparsity Profile", x = "Log10(Mean abundance + 10^-5)", y = "Fraction of zero values - sparsity")
ggsave("plots/p17.png")
# It displays a more continuous and gradual decline in sparse features within the data, but still most sparseness is located towards the left of low abundance


# Another plot to look into especially heteroscedasticity - phenomenon where the variance of microbial abundance is not constant, but changes systematically
# Very often microbial data is right skewed, sparse, and heavily dispersed, which is sort of seen in the previous abundance vs zero fraction plot.

p18 <- ggplot(qc_metrics, aes(x = Log_Mean, y = log10(Variance + 1))) + geom_point(alpha = 0.4, colour = "purple") + labs(title = "Mean and Variance relationship - heteroscedasticity", x = "Log 10(Mean abundance + 1)", y = "Log 10(Variance + 1)")
ggsave("plots/p18.png")
# We see that the variance is increasing with mean abundance, so a positively increasing relationship indicating the heteroscedasticity, there is much more variance with more highly abundant species
p19 <- ggplot(qc_metrics, aes(x = Log_Mean_p, y = log10(Variance + 10^-5))) + geom_point(alpha = 0.4, colour = "purple") + labs(title = "Mean and Variance relationship - heteroscedasticity", x = "Log 10(Mean abundance + 10^-5)", y = "Log 10(Variance + 10^-5)")
ggsave("plots/p19.png")
# The same plot just with a smaller pseudo-count for higher biological accuracy
# The positively increasing relationship is much clearer now, this cone shape of increasing variance with abundance

# With some reading, since we have non-constant variance, more of a spread, statistical method assumptions break, including ANOVA 

# Additionally, since we have removed ~ 80% from previous filtering, sparse filtering may not be a good idea because the previous plots show that a lot of data still contains 0 inflated values;
# Therefore, we could risk removing biologically relevant either low abundant species or true biologically absent samples. 

# With a deeper web search, rCLR - robust Centered Log Ratio transformation could be better suited for my highly sparse data.
# We can experiment to see which gives better results to support downstream analysis;
# While CLR most likely will fail unless imputation is completed because it cannot perform log(0) - I will stick to rCLR;

if(!require("compositions")) install.packages("compositions") # necessary package to complete CLR transformation
library(compositions)
packageVersion("compositions") # Version: 2.0.9

seq_asv_f_clr <- clr(seq_asv_f_mat) # CLR transforming my prevalence and abundance filtered ASV data;

# sanity check;
dim(seq_asv_f_mat)
dim(seq_asv_f_clr)
# They are the same dimension hence all information has been contained

sum(seq_asv_f_clr == 0) # 164437 entries are still 0, hence this would call for some imputation method
sum(is.na(seq_asv_f_clr)) # nicely, there are no NAs in the data

# Based on these small sanity checks, it is preferable to utilise rCLR in this context, why? => rather than imputing data and possibly risking potential artificial data spikes and data distortion, rCLR can handle 0 values more robustly.

if(!require("vegan")) install.packages("vegan")
library(vegan)
packageVersion("vegan") # Version: 2.7.6

# completes rCLR transformation through the vegan package
seq_asv_f_rclr <- decostand(seq_asv_f_mat, method = "rclr")

# sanity check;
dim(seq_asv_f_mat)
dim(seq_asv_f_rclr)
# All clear

# Additionally, one can check the sparsity couunt in the data; this will reveal if the zero counts have been tended to;
sum(seq_asv_f_rclr == 0)
# The result should be and is 0, hence the sparsity has been accounted for;

# HOW? rCLR works as CLR but leaves out any 0 values in the geometric mean, labeling them as NAs. Since PCA and other methods fails with NA - a low rank singular value decomposition is applied, effectively reconstructs a dense matrix without introducing pseudo-count bias.
sum(is.na(seq_asv_f_rclr)) # This shows that there is no NA in the data and will be okay to continue to downstream methods such as PCoA and PCA.

# Visualising the newly transformed and filtered data is the next step to see any underlying behaviour in the data that may still persists;

#------------
# CLR vs rCLR
#------------

# Combine the data to have a consistent data frame to plot;
# Firstly, transform the data horizontally (flip) to receive a long format of the data

if(!require("tidyverse")) install.packages("tidyverse") # package required for the horizontal transformation: pivot_longer()
library(tidyverse)
packageVersion("tidyverse") # Version: 2.0.0

clr_long <- as.data.frame(seq_asv_f_clr) %>%
  pivot_longer(cols = everything(), names_to = "ASV", values_to = "Value") %>% # label the columns for the new data frame as you want
  mutate(Transformation = "CLR")

rclr_long <- as.data.frame(seq_asv_f_rclr) %>%
  pivot_longer(cols = everything(), names_to = "ASV", values_to = "Value") %>%
  mutate(Transformation = "rCLR")

clr_rclr_df <- bind_rows(clr_long, rclr_long) %>%
  filter(!is.na(Value)) # remove any NA values if there appear any, just for the plotting

p20 <- ggplot(clr_rclr_df, aes(x = Transformation, y = Value, fill = Transformation)) + geom_violin(alpha = 0.4) + labs(title = "Difference in distribution W.R.T the transformation (CLR vs rCLR)", x = "Transformation", y = "ASV Values") + facet_wrap(~ Transformation, scales = "free_x") + theme(legend.position = "none") # W.R.T - With Respect To | "free_x" removes the x axis label - zooms into the plots better
ggsave("plots/p20.png")
# This shows the difference in the data distribution after transforming the data. We can see that rCLR has a more symmetrical data spread and not such a dense concentration around the origin as CLR has.
# Additionally, CLR creates the data to be much more right skewed, which is not so adaptive to many downstream statistical testing, including ANOVA.

# PCoA will the choice of visual display - it does not rely on standard euclidean distance matrix, and is considered a gold standard for beta-diversity in microbiome data visualisation
# The distance metric to explore the data with will be "Bray-Curtis" and possible to explore URF - Unsupervised Random Forest to test how close the data points are to each other for clustering charactersitics;


#-------------------------------------------------
# PCoA - Principle Coordinate Analysis (Aitchison)
#-------------------------------------------------

# Since alpha-diversity has been visualised before with library depth, sequencing, and abundance, we can use PCoA to look into the beta-diversity from the data in different individuals;
# A small sanity check to see if the samples have been preserved;
table(colnames(seq_asv_f_mat) == colnames(seq_asv_f_clr))
# outputs the same dimension of features as we expected (474) and all values are true so the samples are still preserved!

# With a further review, Aitchison distance could also favor the data instead of Bray-curtis
# Bray-Curtis is supposed to operate on relative proportions and is vulnerable to compositional bias, while Aitchison respects compositional constraints.
# While I believe the filtering and rCLR transformations will have removed some compositional bias and is it not entirely clear which will perform better;

aitchison_distance <- dist(seq_asv_f_rclr, method = "euclidean") # compute the Aitchison distance between points in the euclidean space
aitchison_pcoa <- cmdscale(aitchison_distance, k = 2, eig = TRUE) # k = 2 is the dimension of the dataset, eigenvalues are true supporting scree plots

# Since the PCoA metrics have been computed, we can construct a data frame to contain our data to them plot it;
aitchison_df <- data.frame(PCo1 = aitchison_pcoa$points[, 1], PCo2 = aitchison_pcoa$points[, 2], Treatment = seq_meta$Treatment, Subject_ID = seq_meta$ID) %>% # collecting all points for PC1 from first column and PC2 from second column respectively
  mutate(distance_centre = sqrt((PCo1 - mean(PCo1))^2 + (PCo2 - mean(PCo2))^2), is_outlier = distance_centre > (mean(distance_centre) + 3 * sd(distance_centre))) # calculate the distance from the center of the data to find the majority cluster and set a threshold of 3 standard deviations to be an outlier for the proposed data

# distance_centre is the euclidean distance from the center of the centroid to the data points, following the Pythagorean theorem
# This helps us define the distance of points and label outliers with distance from the centre of mean with standard deviations.

# We can also add the variance of the data explained to the plots to see the contribution of PCs;

aitchison_var_PCo1 <- round((aitchison_pcoa$eig[1] / sum(aitchison_pcoa$eig[aitchison_pcoa$eig > 0])) * 100, 1)
aitchison_var_PCo2 <- round((aitchison_pcoa$eig[2] / sum(aitchison_pcoa$eig[aitchison_pcoa$eig > 0])) * 100, 1)
# We also now have information that over 80% of the variance is explained already by PC1 and PC2 combined

p21 <- ggplot(aitchison_df, aes(x = PCo1, y = PCo2, colour = Treatment)) + geom_point(aes(shape = is_outlier), size = 1, alpha = 0.8) + stat_ellipse(type = "t", linetype = 2, level = 0.95) + geom_text(data = filter(aitchison_df, is_outlier), aes(label = Subject_ID), vjust = -1, colour = "black", size = 1) + labs(title = "Aitchison distance PCoA (rCLR)", x = paste0("PCo1 (", aitchison_var_PCo1, "%)"), y = paste0("PCo2 (", aitchison_var_PCo2, "%)"))
ggsave("plots/p20.png")
# provides the PCoA plot of PC1 vs PC2 with treatment of diet as the colour coding. The data is mean centered and shows outliers that are 3 sd away from the center, with an ellipse encapsulating the majority of data
# from the plot one can see 2 direct outliers, but we can verify this with the table summary;
table(aitchison_df$is_outlier)
# 2 are returned with the value TRUE, hence we do only have 2 outliers

# Additionally, there is no direct cluster pattern to be seen for the given filtered samples

#-------------------------------------
# Scree plot - PCs variance explained
#-------------------------------------

fviz_pcoa <- as_factoextra_pca(ind.coord = aitchison_pcoa$points, eig = aitchison_pcoa$eig) # convert the pcoa object from cmdscale via this function utilising the eigenvalues = TRUE arguement
p22 <- fviz_eig(fviz_pcoa, addlabels = TRUE, barfill = "blue", barcolor = "blue", linecolor = "red", main = "Scree plot showcasing the PCoAs variance explained (Aitchison)", ylim = c(0, 46))
ggsave("plots/p21.png")
# This proves and displays that all the PCs variance (100% )explained with the first 3 PCoA1, 2, and 3;

# We can also see the connection between other PCs;
aitchison_PCo3 <- cmdscale(aitchison_distance, k = 3, eig = TRUE)
# PC1 vs PC3
aitchison_df_PCo3 <- data.frame(PCo1 = aitchison_PCo3$points[, 1], PCo3 = aitchison_PCo3$points[, 3], Treatment = seq_meta$Treatment, Subject_ID = seq_meta$ID) %>% # collecting all points for PC1 from first column and PC2 from second column respectively
  mutate(distance_centre = sqrt((PCo1 - mean(PCo1))^2 + (PCo3 - mean(PCo3))^2), is_outlier = distance_centre > (mean(distance_centre) + 3 * sd(distance_centre))) # calculate the distance from the center of the data to find the majority cluster and set a threshold of 3 standard deviations to be an outlier for the proposed data

# distance_centre is the euclidean distance from the center of the centroid to the data points, following the Pythagorean theorem
# This helps us define the distance of points and label outliers with distance from the centre of mean with standard deviations.

# We can also add the variance of the data explained to the plots to see the contribution of PCs;

aitchison_var_PCo3 <- round((aitchison_PCo3$eig[3] / sum(aitchison_PCo3$eig[aitchison_PCo3$eig > 0])) * 100, 1)
# We also now have information that over 80% of the variance is explained already by PC1 and PC2 combined, which is already very strong values

p23 <- ggplot(aitchison_df_PCo3, aes(x = PCo1, y = PCo3, colour = Treatment)) + geom_point(aes(shape = is_outlier), size = 1, alpha = 0.8) + stat_ellipse(type = "t", linetype = 2, level = 0.95) + geom_text(data = filter(aitchison_df_PCo3, is_outlier), aes(label = Subject_ID), vjust = -1, colour = "black", size = 1) + labs(title = "Aitchison distance PCoA (rCLR)", x = paste0("PCo1 (", aitchison_var_PCo1, "%)"), y = paste0("PCo3 (", aitchison_var_PCo3, "%)"))
ggsave("plots/p22.png")
# No distinct cluster either, more outliers, and much less explained variance. So from my judgement I will stick to PC1 and PC2 for the majority for explained variance of my data

#======
bray_curtis_distance <- vegdist(seq_asv_f_rclr, method = "bray") # vegdist = dist, just from the vegan package, that contains the Bray-Curtis distance metric
# Running this gives a warning that "results may be meaningless because data have negative entries" 
# This is true after rCLR, I have got rid of any sparsity in my data, but the log ration transformation produced negatives results and therefore Bray-Curtis cannot be used.
#======

# If I want to test Bray-Curtis, then I have to go back and choose a different transformation for the filtered data i.e. TSS - Total Sum Scaling
# https://pmc.ncbi.nlm.nih.gov/articles/PMC9238373/#s2 this paper used both Aitchison and Bray-Curtis, based on Figure 3 on the paper overall Aitchison seems to perform better than Bray-Curtis.
# Additionally, another to check and experiment can be Jaccard distance measure, since it also have a higher AUC-ROC for performance compared to both the previous metrics.


#-----------------------------------------------
# PCoA - Principle Coordinate Analysis (Jaccard)
#-----------------------------------------------

# NOTE: the Jaccard distance metric works with binary numbers, hence we need to convert our abundances to 0 and 1 = non-abundant & abundant
seq_asv_binary <- ifelse(seq_asv_f_rclr > 0, 1, 0)
# Similarly, calculate the Jaccard distance as before; (with vegdist - from the vegan package)
jaccard_distance <- vegdist(seq_asv_binary, method = "jaccard")
jaccard_pcoa <- cmdscale(jaccard_distance, k = 2, eig = TRUE, add = TRUE) # add = TRUE, handles non-euclidean negative eigenvalues, can be seen more in 'Help' in RStudio

jaccard_df <- data.frame(PCo1 = jaccard_pcoa$points[, 1], PCo2 = jaccard_pcoa$points[, 2], Treatment = seq_meta$Treatment, Week = as.factor(seq_meta$time), Subject_ID = seq_meta$ID) # gather the information from score points and meta data to construct a jaccard data frame viable to use for plotting

p24 <- ggplot(jaccard_df, aes(x = PCo1, y = PCo2, colour = Treatment, shape = Week)) + geom_point(size = 1, alpha = 0.8) + stat_ellipse(aes(group = Treatment), type = "t", linetype = 2) + labs(title = "Beta diversity Jaccard distance PCoA (rCLR)", subtitle = "Exploring microbial community shifts across treatment and time", x = "PC1", y = "PC2")
ggsave("plots/p23.png")
# THIS SHOWS THAT JACCARD INDEX CANNOT BE PLOTTED ACROSS rCLR data because we require binary values, and negative values from rCLR will distort our analysis as seen above;

seq_asv_binary_raw <- ifelse(seq_asv_f_mat > 0, 1, 0) # USE THE "RAW" MATRIX SEQ DATA
# Similarly, calculate the Jaccard distance as before; (with vegdist - from the vegan package)
jaccard_distance_raw <- vegdist(seq_asv_binary_raw, method = "jaccard")
jaccard_pcoa_raw <- cmdscale(jaccard_distance_raw, k = 2, eig = TRUE, add = TRUE) # add = TRUE, handles non-euclidean negative eigenvalues, can be seen more in 'Help' in RStudio

jaccard_df_raw <- data.frame(PCo1 = jaccard_pcoa_raw$points[, 1], PCo2 = jaccard_pcoa_raw$points[, 2], Treatment = seq_meta$Treatment, Week = as.factor(seq_meta$time), Subject_ID = seq_meta$ID) %>% # gather the information from score points and meta data to construct a jaccard data frame viable to use for plotting
  mutate(distance_centre = sqrt((PCo1 - mean(PCo1))^2 + (PCo2 - mean(PCo2))^2), is_outlier = distance_centre > (mean(distance_centre) + 3 * sd(distance_centre))) # find potential outliers
  
jaccard_var_PCo1 <- round((jaccard_pcoa_raw$eig[1] / sum(jaccard_pcoa_raw$eig[jaccard_pcoa_raw$eig > 0])) * 100, 1)
jaccard_var_PCo2 <- round((jaccard_pcoa_raw$eig[2] / sum(jaccard_pcoa_raw$eig[jaccard_pcoa_raw$eig > 0])) * 100, 1)

p25 <- ggplot(jaccard_df_raw, aes(x = PCo1, y = PCo2, colour = Treatment)) + geom_point(aes(shape = is_outlier), size = 1, alpha = 0.8) + stat_ellipse(aes(group = Treatment), type = "t", linetype = 2) + labs(title = "Beta diversity Jaccard distance PCoA", subtitle = "Exploring microbial community shifts across treatment", x = paste0("PCo1 (", jaccard_var_PCo1, "%)"), y = paste0("PCo2 (", jaccard_var_PCo2, "%)"))
ggsave("plots/p24.png")
# Now this seems much more natural, with some slightly distinctive clusters appearing within the data

# NOTICEABLY, very little explained variance is displayed through PC1 and PC2

# CHECK OUTLIERS;
table(jaccard_df_raw$is_outlier)
# based on the same conditions there are no outliers, so we will stick to remove the 2 outliers from the Aitchison distance metric.

# There is visual overlap of treatment groups and the time stamps. Meaning, there is not a massive difference for the treatment in this study, there is nothing that displays huge separation.
# Additionally, this can be explained through each individual having somewhat a different microbial community baseline - shown through random distribution and dispersion, instead of clear clusters/separation.

# Therefore, because there are no crazy geometric artifacts seen in the plots, I would say it safe to assume there is no batch effect. So no batch correction will be done.
# From both distance measures, no artificial artifacts could be visually seen...

#---------------------------------------
# Jaccard Scree plot - PCs contribution
#---------------------------------------

fviz_pcoa_jac <- as_factoextra_pca(ind.coord = jaccard_pcoa_raw$points, eig = jaccard_pcoa_raw$eig) # convert the pcoa object from cmdscale via this function utilising the eigenvalues = TRUE arguement
p26 <- fviz_eig(fviz_pcoa_jac, addlabels = TRUE, barfill = "blue", barcolor = "blue", linecolor = "red", main = "Scree plot showcasing the PCoAs variance explained (Jaccard)", ylim = c(0, 10))
ggsave("plots/p25.png")
# This proves and displays that all the PCs variance (100% )explained with the first 3 PCoA1, 2, and 3;


#===================
# Removing Outliers
#===================

# identify the rows for sample outliers (not necessary, but to see their own mini dataset);
outlier_4031 <- seq_asv_f_rclr[seq_meta$ID == "4031", ]
outlier_4018 <- seq_asv_f_rclr[seq_meta$ID == "4018", ]

# combine it and remove it;
seq_asv_f_cleaned <- seq_asv_f_rclr[!(seq_meta$ID %in% c("4031", "4018")), ] # dropping the outlier samples with their duplicates
# sanity check to see if it was successful, since are removing 2 samples but 6 time points for each sample, we expect "12 samples" to be removed;
dim(seq_asv_f_rclr) - dim(seq_asv_f_cleaned) # 484 - 472 = 12

# I will also stick with only removing the 2 samples above, not to weaken downstream analysis with a large reduction of sample numbers


#===================
# PERMAVONA Testing
#===================


# https://www.zubairkhalid.com/knowledge/bioinformatics/permanova-in-microbiome-studies-assumptions-implementation-and-reporting-standards
# We want to test the differences in microbial community composition between the two groups, done by partitioning variation in a distance matrix.
# PERMANOVA is sensitive to both differences in group location - centroids, and differences in group dispersion - spread.
# Due to this sensitivity, one should check for homogeneity of multivariate dispersion.

# Overall we want to test here if FV as treatment induces a significant shift in the gut microbiome when compared to UD over time.
# NOTE: since we have "repeat" individual samples - can use restricted permutations, so that samples are shuffled within each individual preserving the baseline microbiome
# NOTE: homogeneity - to test whether group centroids are different, if clustering of one group and scattering of another occurs then the test inflates results... > test with PERMDISP2

#-------------
# Homogeneity
#-------------

# calculate the average distance of sample points to their group median, to notify if the dispersion is significant.
# p > 0.05 = assumption of homogeneity is met, P < 0.05 = assumption is not met and it is significant

# Make all meta data variables factors, else error;
seq_meta$Treatment <- as.factor(seq_meta$Treatment)
seq_meta$time <- as.factor(seq_meta$time)
seq_meta$ID <- as.factor(seq_meta$ID)

# Jaccard will have to be used here for presence-absence patterns - therefore, we cannot use rCLR data we have before, but the seq_asv_filtered, since I came to the conclusion that sparsity is not a damaging factor, since most was already filtered out.

# Before, we must define a non-rCLR transformation and add an ID column to understand our sample names;
remove_outliers <- c("4031", "4018")
keep_IDs <- which(!(seq_meta$ID %in% remove_outliers))

seq_asv_f_cleaned_no_rclr <- seq_asv_f_mat[keep_IDs, ]
seq_meta_cleaned <- seq_meta[keep_IDs, ]

# sanity check;
dim(seq_asv_f_mat) - dim(seq_asv_f_cleaned_no_rclr) # 12 again - worked
dim(seq_meta) - dim(seq_meta_cleaned) # also 12

# We have already Jaccard distance, just with the cleaned version (no rCLR) for rigorousness and to avoid downstream conflicts completed again;
seq_asv_bin <- ifelse(seq_asv_f_cleaned_no_rclr > 0, 1, 0) # USE THE "RAW" MATRIX SEQ DATA THAT IS REMOVED OF OUTLIERS
# Similarly, calculate the Jaccard distance as before; (with vegdist - from the vegan package)
jaccard_dist <- vegdist(seq_asv_bin, method = "jaccard")
# sanity check;
attr(jaccard_dist, "Size")
nrow(seq_meta_cleaned)
# dimensions must match after cleaning and filtering, and they do

# We have our data prepared for PERMANOVA and sub-set the correct samples for processing - "cleaned" data frames

treatment_dispersion <- betadisper(jaccard_dist, seq_meta_cleaned$Treatment) # betadisper - multivariate homogeneity of groups dispersions test
homogeneity_anova <- anova(treatment_dispersion)
print(homogeneity_anova)
# PERMUTATIONS are more robust than a standard ANOVA;
homogeneity_perm <- permutest(treatment_dispersion, permutations = 999) # is more robust than ANOVA, due to the complexity of microbiome studies and the violation of assumptions, permutest is a better option
print(homogeneity_perm) # even though the p value is higher, it tests a permutation distribution of the F-statistic under the null-hypothesis
homogeneity_perm_df <- data.frame(homogeneity_perm$tab) # just for an easier visual

# homogeneity - Pr value = 0.00003725 < 0.05 so it is statistically significant, meaning we can accept that there is a significant difference comparing the two groups Control vs Treated for the spread (variance in data points).
# We need to now understand where this difference of dispersion comes from;

# Since we have significant dispersion p value results, we must tread carefully, we have more of a heterogeneous spread - so the results could be driven by some unknown spread, and not necessarily the true shift we are looking for
# Both tests have a significant p value and hence, the treatment and control groups have different levels of dispersion.

p27 <- boxplot(treatment_dispersion, main = "Dispersion Spread UD vs FV", xlab = "Group", ylab = "Distance to centroid")
ggsave("plots/p26.png")
# Although, now this plot contradicts what is mentioned above, since the difference spread is visually non-existent between the two groups.
# This can be due to the high degrees of freedom regarding the samples, that any "microscopic" difference can be detected as significant.

#----------------------
# Restricted PERMANOVA
#----------------------

# We want to use restricted to only shuffle individual's samples via how()
# This will allow us to test whether the centroids are distinct and if they are pushed apart in space
# Restricted permanova, restricting on unique subject-block permutations to preserve the individual sample integrity - 999 permutations, chosen randomly but enough times.
PERMANOVA_control <- how(block = seq_meta_cleaned$ID, nperm = 999)

# Analysis of variance using distance matrices - looking for the underlying Treatment and Week relationship via multiplicative equation
PERMANOVA_model <- adonis2(jaccard_dist ~ Treatment * time, data = seq_meta_cleaned, permutations = PERMANOVA_control)
print(PERMANOVA_model)

# Outputs a p value of 0.001 => significant results, since we looked at the interaction of treatment and time, we want to know if based on the treatment type the difference in microbial community is significant or not.
# Here, evidently it shows it is, if we take p <= 0.05 as a significance threshold. We deduce that FV successfully altered the individuals microbial composition over time compared to the UD.

# NOTE: since we have been using Jaccard distance measures - it is a global view of the data, we cannot see the macro-level community shifts...
# To see which specific species of bacteria could be potentially initiating certain shifts.
# ANCOM-BC2 - Analysis of Compositions of Microbiomes with Bias Correction 2 -> can identify individual ASVs that change significantly over time and also correcting for any sequencing bias that may be present.

# Moreover, ANCOM-BC2 required raw, filtered, and un-transformed data -> no rCLR! Similarly as for Jaccard. Else it cannot correct for sequencing bias of course.

# BiocManager::install("ANCOMBC") # package and library for ANCOM-BC2
# library(ANCOMBC)
# packageVersion("ANCOMBC") # Version: 2.14.0

# https://bioconductor.posit.co/packages/3.19/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html
# Information regarding the use of ANCOMBC function;

# Allows one to see which species of the microbiome absolute abundance, per unit volume of the whole composititon are significantly different with changes in the covariate of interest (groups - UD/FV)
# An Error is given if seq_asv_f_cleaned_no_rclr is not applied with as.matrix(), even though it is a numeric data frame

## seq_asv_count_mat <- as.matrix(seq_asv_f_cleaned_no_rclr)
## storage.mode(seq_asv_count_mat) <- "numeric"
## ANCOMBC_output <- ancombc2(data = seq_asv_count_mat, taxa_are_rows = FALSE, assay_name = "counts", meta_data = seq_meta_cleaned, tax_level = NULL, fix_formula = "Treatment * time", rand_formula = "(1 | ID)", p_adj_method = "holm", lib_cut = 0, group = "Treatment", struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, n_cl = 1)

# 1) select data to perform ancombc2 on the matrix with no rCLR
# 2) taxa are columns - the species
# 3) select the count table in the data
# 4) select the cleaned meta data
# 5) We do not have species hierarchy - levels i.e. genus
# 6) same behaviour formula as before to test treatment and time
# 7) random chosen intercept formula for repeated measures in the individuals
# 8) multiple testing correction method - holm, others exist such as, Bonferroni, FDR, etc
# 9) set to 0 because we already completed filtering pre-hand
# 10) the main grouping variable, it is what we want to test
# 11) since the data is not transformed, we want the function to detect 0s
# 12) conservative variance estimation
# 13) alpha error, significance threshold
# 14) number of CPU cores to use, this device has 8 cores, and the dataset is not large after filtering. So I propose to just use 1.

# I am running into a data frame error, that seems for now un-fixable, so I will skip it for now.
# I believe my sample data is consistent between the data frames and any previous fixes did not change the error message or improve anything...

# What would of been nice to extract from this multivariate test, is to see which species of the microbiome community are responding directly to the treatment effect of FV over the period of 8 weeks.


#==================================================================
# RM-ASCA - Repeated Measures ANOVA Simultaneous Component Analysis
#==================================================================


# # For RM-ASCA we will use our filtered, rCLR transformed data, and cleaned data frame;
# 
# # Sanity check to see the rCLR data distribution;
# rclr_long <- as.data.frame(seq_asv_f_cleaned) |> tidyr::pivot_longer(everything(), names_to = "ASV", values_to = "rCLR")
# ggplot(rclr_long, aes(x = rCLR)) + geom_histogram(bins = 50) + theme_classic() + labs(title = "ASV distribution after rCLR")
# # This definitely shows a normal distribution, but a very high concentration close to the origin
# 
# # Also, to create our design matrix - this incorporates the individual's ID, time, and treatment together for further analysis in a linear model; (y = mx + c)
# ASCA_df <- as.data.frame(seq_asv_f_cleaned) # data to use for ASCA, cleaned and rCLR
# ASCA_design <- model.matrix(~ 0 + seq_meta_cleaned$Treatment * seq_meta_cleaned$time + seq_meta_cleaned$ID)
# # This will allow the mathematical separation of of systematic experimental variations from individual variations.
# # We want to capture that potential heterogeneous individual response
# 
# # First decompose the fitted effects per factor and extract the fitted values for the equation Treatment * time + ID interaction;
# ASCA_fit_model <- lm(as.matrix(ASCA_df) ~ ASCA_design)
# 
# coef_mat <- coef(ASCA_fit_model)
# interaction_rows <- grep("Treatment.*time", rownames(coef_mat))
# 
# # Extract the residuals and the effect matrices 
# ASCA_effect_matrix <- coef_mat[interaction_rows, ]
# 
# # For visualisation, I will lean for PCA to visualise underlying trajectories of the FV treatment from the interaction effect matrix;
# ASCA_interaction_pca <- prcomp(ASCA_effect_matrix, scale. = FALSE)
# # Create the necessary data frame for the trajectory;
# ASCA_plot_df <- data.frame(PC1 = ASCA_interaction_pca$x[, 1], PC2 = ASCA_interaction_pca$x[, 2], Treatment = seq_meta_cleaned$Treatment, Time = seq_meta_cleaned$time, Subject_ID = seq_meta_cleaned$ID)
# # PC1 = 1st column of the PC points, PC2 the 2nd column
# 
# # WARNING: this only produced 1 column that can be extracted as the 1st PC, we need a different solution
# 
# ggplot(ASCA_plot_df, aes(x = PC1, y = PC2, colour = Treatment, shape = Time)) + geom_point(size = 1, alpha = 0.8) + geom_line(aes(group = Subject_ID), alpha = 0.3, colour = "grey50") + stat_ellipse(aes(group = Treatment), type = "t", linetype = 2) + labs(title = "RM-ASCA Trajectory Analysis", subtitle = "Track individual shifts with respect to time * treatment", x = "ASCA component 1", y = "ASCA component 2")
# 
# # Again, we still see major overlap between the two treatment groups. 

#----------

# Previous, work was not successful due to a data frame being out of bounds, this comes from the prcomp function creating a 1 dimensional data object that does not have a 2nd column for a 2nd PC
# Different solution;
ASCA_df <- ASCA_df <- as.data.frame(seq_asv_f_cleaned)
ASCA_design <- model.matrix(~ seq_meta_cleaned$Treatment * seq_meta_cleaned$time) # create a design matrix with an interaction equation looking at relationship between the treatment and time
ASCA_fit_model <- lm(as.matrix(ASCA_df) ~ ASCA_design) # fitting a linear regression model on the cleaned ASV matrix accounting for our design matrix interaction equation
interaction_cols <- grep("Treatment.*time", colnames(ASCA_design)) # collect the interaction column names

interaction_component_matrix <- ASCA_design[, interaction_cols, drop = FALSE] %*%
  coef(ASCA_fit_model)[grep("Treatment.*time", rownames(coef(ASCA_fit_model))), , drop = FALSE]

ASCA_int_pca <- prcomp(interaction_component_matrix, scale. = FALSE) # Principal component extraction from the interaction componenet matrix

ASCA_plot_df <- data.frame(PC1 = ASCA_int_pca$x[, 1], PC2 = if(ncol(ASCA_int_pca$x) > 1) ASCA_int_pca$x[, 2] else 0, Treatment = seq_meta_cleaned$Treatment, Time = seq_meta_cleaned$time, Subject_ID = seq_meta_cleaned$ID)
# For PC2 in case there are values missing or the column does not exist it will choose the option that should work

p28 <- ggplot(ASCA_plot_df, aes(x = PC1, y = PC2, colour = Treatment, shape = Time)) + geom_point(size = 2, alpha = 0.8) + geom_line(aes(group = Subject_ID), alpha = 0.3, colour = "grey50") + stat_ellipse(aes(group = Treatment), type = "t", linetype = 2) + labs(title = "RM-ASCA Trajection analysis", subtitle = "Track individual shifts with repsect to interaction components", x = "ASCA component 1", y = "ASCA component 2")
ggsave("plots/p27.png")
# The second component is so infinitely small that the interaction behaviour cannot be seen. This can be due to the 1 degree of freedom given by time, and it seems that the second component has practically no explained variance.
# Therefore, it could be recommended to plot on a 1 dimensional scale i.e. Time vs PC1;

ASCA_plot_df$Time <- as.factor(ASCA_plot_df$Time)
ASCA_plot_df$Treatment <- as.factor(ASCA_plot_df$Treatment)

# Provide individual trajectories with geom_line and individual data points on the plot with geom_point;
# Also group mean trend and group mean points with stat_summary respectively;
p29 <- ggplot(ASCA_plot_df, aes(x = Time, y = PC1, colour = Treatment, group = Subject_ID)) + geom_point(aes(shape = Time), size = 2) + geom_line(alpha = 0.3, linewidth = 0.8) + stat_summary(aes(group = Treatment), fun = mean, geom = "line", linewidth = 1.5, linetype = 1) + stat_summary(aes(group = Treatment), fun = mean, geom = "point", size = 4) + labs(title = "RM-ASCA Interaction Trajectory Profile", subtitle = "Individual subject paths and group mean trends along Component 1", x = "Time Point", y = "Interaction Component 1 Score") + 
# The plot itself doesn't provide much information
# Most likely the data was compressed and not expressed as I would of liked to
# But at the end of the day this plot highlights that there actually is some sort of distinctive temporal shift within the 8 weeks for the FV group compared to the control UD group

# Maybe we can adjust the plot to make it look nicer or more concrete;
  facet_wrap(~ Treatment, scales = "free_y") + theme_minimal(base_size = 12)
ggsave("plots/p28.png")

#========================================================================
# Machine Learning - Are there distinct patterns in microbiome responses
#========================================================================


# Collect all the data that was originally split up before to have 1 consistent table;
# NOTE: Since inflammatory markers only have unique ID numbers and therefore less row numbers, we need to complete a left join to accept the data available correctly
meta_df <- seq_meta_cleaned %>%
  mutate(ID = as.factor(ID))

inf_markers$ID <- as.factor(inf_markers$ID) # cannot join factors and integers, so must convert prior
meta_inf_combined <- meta_df %>%
  left_join(inf_markers, by = "ID")

asv_df <- as.data.frame(seq_asv_f_cleaned) # for consistency

ml_complete_data <- cbind(meta_inf_combined, asv_df) # collect together all necessary data to confirm everything is in place
# sanity check;
dim(ml_complete_data) # now the dimensions are consistent with our upstream analysis

# Deciding here that I will not need the columns 4:18, to make the machine learning easier by not providing unnecessary data
drop_columns <- c("observed_features", "chao1", "simpson", "Baxis1", "Baxis2", "Baxis3", "Jaxis1", "Jaxis2", "Jaxis3", "Waxis1", "Waxis2", "Waxis3", "Uaxis1", "Uaxis2", "Uaxis3")
machine_data <- ml_complete_data %>%
  select(-one_of(drop_columns)) %>%
  mutate(ID = as.factor(ml_complete_data$ID), time = as.factor(ml_complete_data$time), Treatment = as.factor(ml_complete_data$Treatment))
# Selecting all the variables and data frames, now we have all the correct format and data for the next steps

# Since we want to predict the treatment group, we will set our target variable as the treatment group;
  
# Sanity check for NAs, since ML breaks with missing values;
which(is.na(machine_data)) # outputs where the NAs lie
# They seem to be only in the inflammatory markers data frame, we cannot just remove them without ruining the original data composition and cannot just impute 0 instead of NA
# A solution could be to complete MAR - Missing At Random imputation, as I cannot would assume from the blood tests the inflammation was not picked up and telling a machine learning algorithm that something is missing or 0 will allow it to be mislead.
# It is not mentioned in the paper how the values became missing - but likely that MAR applies here

# Since we don't want to remove samples, lose dimensions and statistical power, we must impute the values

# To test machine learning algorithms, XGBoost, Elastic Net and Random Forest are familiar
# NOTE: XGBoost can handle NAs, while the other two cannot
# Therefore, for a strong imputation algorithm, I will use MissForest - an iterative Random Forest imputation algorithm , handling non-linear relationships without assuming any specific data distributions.

# Necessary for machine learning data splitting
if(!require("caret")) install.packages("caret")
library(caret)
packageVersion("caret") # Version: 7.0.1


#===================================================

# # I realised I completed this without specifying the unique samples
# train_id <- createDataPartition(machine_data$Treatment, p = 0.8, list = FALSE) # partition the data into a 80% and 20% split
# train_raw <- machine_data[train_id, ] # split into training and testing data that the algorithms use to predict and then validate on respectively;
# test_raw <- machine_data[-train_id, ]
# 
# # for missForest, the factor variables need to be removed, only the numerical data inflammatory and ASVs are required
# drop_columns_1 <- c("ID", "Treatment")
# 
# # subsetting the train and test data for only the numerical values
# X_train_miss_forest <- train_raw[, !colnames(train_raw) %in% drop_columns_1]
# X_test_miss_forest <- test_raw[, !colnames(test_raw) %in% drop_columns_1]
# 
if(!require("missForestPredict")) install.packages("missForestPredict") # package for missForest
library(missForestPredict)
packageVersion("missForestPredict") # Version: 1.0.1
# 
# machine_imputation <- missForest(xmis = X_train_miss_forest, maxiter = 10, ntree = 100)
# # x_mis = the data containing missing values
# # maxiter = maximum number of iterations, but the algorithm will stop automatically if it converges
# # ntree = number of trees
# 
# # Ran way over 30 mins....
# 
# X_train_clean <- machine_imputation$ximp # extract the imputed training figures only
# 
# # Now we need to validate the results;
# # To test the trained data now onto the test data - with no leakage;
# X_test_clean <- missForestPredict(machine_imputation, newdata = X_test_miss_forest)
# 
# # Imputation was completed on the train data, validated on the test data and now ready to be extract as a matrices to be used for ML;
# X_train_final <- as.matrix(X_train_clean) # These are our train-test split data for ASV and inflammation respectively
# X_test_final <- as.matrix(X_test_clean)
# 
# # Additionally, define our target for our predictors to "search" for;
# Y_train_final <- train_raw$Treatment
# Y_test_final <- test_raw$Treatment

#===================================================

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Previous faulty work, but the same pipeline
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


set.seed(1) # set this for exact reproducibility

# Before any machine learning data splitting should almost 100% completed, to remove any process of data leakage
# This can be done with a train-test split. Usually 80-20 is chosen, where the test set is partitioned to be hidden from the algorithm;

# NOTE: One must be extremely careful since our IDs are not unique and we risk more data leakage within the algorithms. 
# Therefore, we must define unique IDs to separate them accordingly and preserve the correct analysis;
unique_ID <- unique(machine_data$ID)

train_subject_id <- sample(unique_ID, size = floor(0.8 * length(unique_ID))) # partition the unique sample IDs as a 80/20 split

train_uni_raw <- machine_data %>%
  filter(ID %in% train_subject_id) # splitting the rows regard the subject it belongs to
test_uni_raw <- machine_data %>% # everything that was not used by train_uni_raw
  filter(!ID %in% train_subject_id)

# sanity check;
nrow(train_uni_raw) + nrow(test_uni_raw) # the total is correct, so this is still good news.

imputation_dropped_cols <- c("ID", "Treatment") # these are columns that are not necessary for the imputation, as it does not help the process and causes further data leakage.
# While we keep "time" here because that plays a big role in understand the pattern that an individual has multiple samples and we want to discoer patterns with time.

X_train_uni_mF <- train_uni_raw[, !colnames(train_uni_raw) %in% imputation_dropped_cols]
X_test_uni_mF <- test_uni_raw[, !colnames(test_uni_raw) %in% imputation_dropped_cols]

# From parallel computing lecture, for imputation taking so long we can make it more efficient, by splitting the task up;
if(!require("doParallel")) install.packages("doParallel")
library(doParallel) # allows the use of parallel programming to split the tasks into smaller chunks and then collect them together for faster operation.
# The reason I include this is because previously the missForest ran around 40 minutes on default settings, and I believe it could be more efficient;
packageVersion("doParallel") # Version: 1.0.17

registerDoParallel(cores = 6) # my system has 8 cores and it is not a good idea to allocate 100% of the cores for this tasks to not run the risk of crashing the system

# Additionally, since I ran missForest at 10 iterations and 100 trees and it was somewhat inefficient with time, I reduce both by 1/2 respectively;
set.seed(1)
# Apply MAR imputation random forest function on training data 
machine_learning_imputation <- missForest(xmis = X_train_uni_mF, maxiter = 5, ntree = 50, parallelize = "forests") # forest calls the specified cores above to be used, in the sense to parallelise the trees created, so the n-trees get split across my 6 cores
# Much faster run time, than with default parameters.

# Validate the imputation and check that the distribution is still normal for the downstream analysis;

#-----------------------------------------
# Plot Imputation - validate distribution
#-----------------------------------------

# First we need the imputed values from the missForest imputed data object;
X_train_imputed <- machine_learning_imputation$ximp

# Pick a column to check the distribution and we can then overlay with the non-imputed data to check the difference
imputed_feature <- colnames(X_train_imputed)[2] # use the first biomarker, chosen arbitrarily

# the dimensions of both the before and after imputed data is the same. Specifically, the length of the vector (check with: length(X_train_uni_mF[[imputed_feature]]) & length(X_train_imputed[[imputed_feature]]))
imputed_density_df <- data.frame(Value = c(X_train_uni_mF[[imputed_feature]], X_train_imputed[[imputed_feature]]), Dataset = rep(c("Before-Imputation (NAs included)", "After-Imputation (via mF)"), c(length(X_train_uni_mF[[imputed_feature]]), length(X_train_uni_mF[[imputed_feature]]))))
# simply just creating a data frame with the train split data used for the random forest imputation and the imputed data, with labels of data with NAs and without before and after imputation respectively
# plot the differences in density;
p30 <- ggplot(imputed_density_df, aes(x = Value, fill = Dataset)) + geom_density(alpha = 0.4) + labs(title = paste("Validating distribution for:", imputed_feature), x = "Feartue Value", y = "Density")
ggsave("plots/p29.png")
#
# NOTE: the density distribution for the before and after data imputation displays near perfect overlapping behaviour meaning that there is no unnatural or artificial spikes and random data imputation for the downstream analysis.
# One can expect more right skewness for inflammatory markers, since the majority of the population have a lower baseline concentration for a biomarker derived from  a circulating gut microbiome (TMAO1 - TriMethylAmine N-Oxide 1 -> functioning as a powerful circulating gut microbiome and active mediator of vascular inflammation)
#

#
# SIDE NOTE: using RF and XGB - right skewness is not a problem, tree based models are NON-PARAMETRIC - models depend more on the median and less on the mean, with no data distribution assumptions
# While, for Elastic Net - a linear model -> right skewness can reduce performance, this can be the reason for the top features extracted from EN are all ASVs and none of the highly skewed inflammatory markers.
#

# plot the density side by side to get a different comparison perspective - separate on the "Dataset; (remove legend for more space)
p31 <- ggplot(imputed_density_df, aes(x = Value, fill = Dataset)) + geom_density(alpha = 0.6) + facet_wrap(~ Dataset, scales = "free_y") + labs(title = paste("Side-by-Side distribution  validation for:", imputed_feature), x = "Feartue Value", y = "Density") + theme(legend.position = "none")
ggsave("plots/p30.png")

#---- data validation completed
# Continue to model training:

X_train_clean <- machine_learning_imputation$ximp # extract the imputed data frame
X_test_clean <- missForestPredict(machine_learning_imputation, newdata = X_test_uni_mF) # validation step on the test data

# Now we can also drop the "time" factor because we only want the algorithms to learn and predict on numerical ASV and inflammatory data

ml_dropped_cols <- c("time")
# Extract the final pieces for training and testing ready for ML with matrix formats
X_train_final <- as.matrix(X_train_clean[, !colnames(X_train_clean) %in% ml_dropped_cols])
X_test_final <- as.matrix(X_test_clean[, !colnames(X_test_clean) %in% ml_dropped_cols])

Y_train_final <- train_uni_raw$Treatment
Y_test_final <- test_uni_raw$Treatment

# Now we have fully filtered, outlier detected, transformed, NA - imputed data completely ready for ML;

# Elastic Net and XGBoost require numeric classes for targets;
y_train_num <- as.numeric(Y_train_final) - 1 # -1 reduced the inputs to a binary input, i.e. in terms of 0 and 1 and not 1 and 2
y_test_num <- as.numeric(Y_test_final) - 1

# We also want to prevent any further cross individual data leakage, so we must set up group folds to preserve the individual ID integrity;
group_folds <- groupKFold(train_uni_raw$ID, k = 5) # it is recommended to split into either group folds of 5 or 10, since splitting into our 86 unique samples will make our mahcine learning algorithms extremely inefficient

# Lastly, we need to convert the folds above into fold IDs so that they can be identified by the algorithms - Elastic Net and XGBoost
# Elastic Net:
fold_id <- rep(0, nrow(X_train_final))
for(i in seq_along(group_folds)) {
  held_out_ids <- setdiff(1:nrow(X_train_final), group_folds[[i]]) # provides a mean of cleanly separating the folds so they become mutually exclusive and cannot be learned from each other
  fold_id[held_out_ids] <- i
}

# groupKFolds gives the separation for TRAINING indices at each fold
# and we identify the HELD-OUT indices at each fold

# XGBoost:
# Converting the caret split data of training folds into validation / HELD-OUT folds;
boost_validation_folds <- lapply(group_folds, function(train_subject_id) {
  setdiff(1:nrow(X_train_final), train_subject_id)
})

#------------------
# EN - Elastic Net
#------------------

if(!require("glmnet")) install.packages("glmnet")
library(glmnet)
packageVersion("glmnet") # Version: 5.1

set.seed(1)

# set results to compare later;
best_AUC <- 0
best_alpha <- 0 
best_EN_fit <- NULL

# We will perform a little grid search also for alpha between 0 and 1, that will loop across and provide the best fit;
for (alpha in seq(0, 1, by = 0.25)) { # loop for alpha from 0 to 1 at 0.25 steps (4 times)
  EN_fit <- cv.glmnet(x = X_train_final, y = y_train_num, family = "binomial", type.measure = "auc", foldid = fold_id, alpha = alpha)
  max_AUC <- max(EN_fit$cvm)
  if (max_AUC > best_AUC) {
    best_AUC <- max_AUC
    best_alpha <- alpha
    best_EN_fit <- EN_fit
  }
}

# Now validate the test data;
# This uses the most optimal lambda - controlling the overall strength of the model penalty (ridge penalty)
EN_prediction <- predict(best_EN_fit, newx = X_test_final, s = "lambda.min", type = "response")[, 1]

#--------------------
# RF - Random Forest
#--------------------

# Here the use of ranger may be better as it is parallel and much faster than the standard random forest;
if(!require("ranger")) install.packages("ranger")
library(ranger)
packageVersion("ranger") # Version: 0.18.0

set.seed(1)
# Train a random forest model, NOTE: permutation importance is ON - this will help discover the key microbial species in the dataset;
RF_model <- ranger(x = X_train_final, y = Y_train_final, num.trees = 500, mtry = floor(sqrt(ncol(X_train_final))), probability = TRUE, importance = "permutation")
# Setting our training data, deciding the number of tree = 500 (default from the help page), mtry - number of variables to split by at each node (formula is default), probability forest, importance - measure for classification, seed for reproducibility

# Validate on the test data and extract the results;
RF_prediction_object <- predict(RF_model, data = X_test_final)
RF_prediction <- RF_prediction_object$predictions[, 2] # extracting the probabilities for the positive class

#---------------
# XGB - XGBoost
#---------------

if(!require("xgboost")) install.packages("xgboost") # package to execute the XGBoost model algorithm
library(xgboost)
packageVersion("xgboost") # Version: 3.2.1.1

set.seed(1)

# This is the necessary matrix data storage for XGBoost;
d_train <- xgb.DMatrix(data = X_train_final, label = y_train_num)
d_test <- xgb.DMatrix(data = X_test_final, label = y_test_num)

# setting XGBoost parameters to work with in the algorithm;
XGB_parameters <- list(objective = "binary:logistic", eval_metric = "auc", max_depth = 4, eta = 0.05, subsample = 0.8, colsample_bytree = 0.8) 
# https://xgboost.readthedocs.io/en/stable/parameter.html

# max_depth - default 6, making it larger may over-fit
# eta - learning rate, default 0.3, smaller step size shrinkage makes the boosting more conservative
# subsample - same as the train/test split 0.8 > subsampling ratio for training instances at each iteration
# col_sample_by_tree - subsample ratio of columns when constructing each tree, same as before.

# We can insert the manual row indices folds to preserve longitudinal independence;
XGB_fit <- xgb.cv(params = XGB_parameters, data = d_train, nrounds = 300, folds = boost_validation_folds, early_stopping_rounds = 20, verbose = FALSE) 
# nrounds - max number of boosting iterations
# early_stopping_rounds - when there is no improvement in the model, this is the amount of boosting rounds that will occur to stop the model.
# verbose - XGBoost will be silent and not output everything in the console.

# doing this prior, saves time for the below real model, since it will find the best XGBoost "fit" - sort of like a grid search

XGB_model <- xgb.train(params = XGB_parameters, data = d_train, nrounds = XGB_fit$early_stop$best_iteration)

XGB_prediction <- predict(XGB_model, newdata = d_test)

#------------------
# Model Evaluation
#------------------

# Using OOS - Out Of Sample testing -> it can evaluate a model on the unseen data and will check how well it generalises, while preventing over-fitting
# Looking at ROC - Receiver Operating Characteristic for our model evaluation;
if(!require("pROC")) install.packages("pROC")
library(pROC)
packageVersion("pROC") # Version: 1.19.1

set.seed(1)
# Retrieve all the ROC for each model on the predicted data extract from each respective model;
ROC_EN <- roc(y_test_num, as.vector(EN_prediction), quiet = TRUE) # to not run into any error downstream, as.vector is usually required
ROC_RF <- roc(y_test_num, RF_prediction, quiet = TRUE)
ROC_XGB <- roc(y_test_num, XGB_prediction, quiet = TRUE)

# Additionally, retrieve all AUC for each respective model;
AUC_EN <- as.numeric(auc(ROC_EN))
AUC_RF <- as.numeric(auc(ROC_RF))
AUC_XGB <- as.numeric(auc(ROC_XGB))

# We have collected our metrics from OOS prediction data, and this will allow us to see how effective the model where on parsing the rCLR data to predict the Treatment Group classification
# Testing different types of models;
# EN - linear and enforces penalties
# RF - non-linear and evaluates importance based on prediction accuracy
# XGB - sequential booster and optimised largely on local feature splits

# Since the data is a combination of ASV and inflammatory markers as the predictors, and the target is a categorical outcome of treatment group, 
# we aim to find how specific this combination relates to an individuals treatment state.
# Moreover, one could potentially deduce also which species and inflammation markers related to either FV and UD for each specific individual - telling a story of the main drivers within a community.


#=========
# Results
#=========

# Evaluating the model's performances;
EN_prediction_mat <- as.factor(ifelse(EN_prediction > 0.5, "FV", "UD"))
RF_prediction_mat <- as.factor(ifelse(RF_prediction > 0.5, "FV", "UD"))
XGB_prediction_mat <- as.factor(ifelse(XGB_prediction > 0.5, "FV", "UD"))

# Build the confusion matrices for each model based on the factors of the treatment groups and so a significant outcome of at least 0.5 of a prediction factor, validated on the target data set
EN_cm <- confusionMatrix(EN_prediction_mat, Y_test_final)
RF_cm <- confusionMatrix(RF_prediction_mat, Y_test_final)
XGB_cm <- confusionMatrix(XGB_prediction_mat, Y_test_final)

# Create a table with all the results;
summary_models <- data.frame(Model = c("Elastic Net", "Random Forest", "XGBoost"), 
                             Accuracy = c(EN_cm$overall["Accuracy"], RF_cm$overall["Accuracy"], XGB_cm$overall["Accuracy"]), 
                             Sensitivity = c(EN_cm$byClass["Sensitivity"], RF_cm$byClass["Sensitivity"], XGB_cm$byClass["Sensitivity"]), 
                             Specificity = c(EN_cm$byClass["Specificity"], RF_cm$byClass["Specificity"], XGB_cm$byClass["Specificity"]), 
                             AUC = c(AUC_EN, AUC_RF, AUC_XGB))
print(summary_models)

# Results are interesting (very weak), XGBoost performed the best according to the AUC, but had the worst accuracy and sensitivity compared to others, while dominating in specificity.
# This means that XGBoost performs the best for this data in specifying the intricate relationships shared with ASV and inflammation markers to treatment groups.

#========

EN_coef <- coef(best_EN_fit, s = EN_fit$lambda.min) # extract the smallest penalty coefficients and create a data frame for EN containing the Feature and its importance;
EN_importance <- data.frame(Feature = rownames(EN_coef)[-1], Importance = abs(as.vector(EN_coef)[-1])) %>%
  arrange(desc(Importance)) %>% 
  head(14)

# We can visualise this data frame to see the resulting effects
p32 <- ggplot(EN_importance, aes(x = reorder(Feature, Importance), y = Importance)) + geom_col(fill = "purple") + coord_flip() + labs(title = "EN Top 14 Features", x = "Feature", y = "|Coefficient|")
ggsave("plots/p31.png")

# Similarly for the Random Forest model;
RF_importance <- data.frame(Feature = names(ranger::importance(RF_model)), Importance = ranger::importance(RF_model)) %>%
  arrange(desc(Importance)) %>%
  head(14)
# I included 14 here as after 12 there where ASVs included in the feature extraction, not just inflammation markers

# FYI: ranger:: is included since there was likely a package conflict

p33 <- ggplot(RF_importance, aes(x = reorder(Feature, Importance), y = Importance)) + geom_col(fill = "red") + coord_flip() + labs(title = "RF Top 14 Features", x = "Features", y = "Importance measure")
ggsave("plots/p32.png")

# Lastly, for XGBoost;
# https://xgboost.readthedocs.io/en/latest/r_docs/R-package/docs/reference/xgb.plot.importance.html
XGB_importance <- xgb.importance(model = XGB_model)
par(mar = c(5, 15, 4, 2)) # margin vectors to adjust the canvas when needed to fit all information
xgb.plot.importance(XGB_importance, top_n = 14, main = "XGB Top 14 Features", xlab = "Importance measure", cex = 0.7, cex.axis = 0.7, left_margin = 17)
title(ylab = "Features", line = 10)

# Choose between the plots as there are the same just different work-around. In my opinion, below is simpler

if(!require("Ckmeans.1d.dp")) install.packages("Ckmeans.1d.dp") # apparently this package is required to plot with ggplot below?
library(Ckmeans.1d.dp)
packageVersion("Ckmeans.1d.dp") # Version: 4.3.6 
p34 <- xgb.ggplot.importance(XGB_importance, top_n = 14) + labs(title = "XGB Top 14 Features", y = "Importance measure", x = "Features") + guides(fill = FALSE) + scale_fill_manual(values = c("blue", "blue", "blue", "blue", "blue"))
ggsave("plots/p33.png")
 # guides removes the legend, so there is more space on the plot itself.
# scale_fill_manual to adjust the colours since it automatically clusters colours - in a non-trivial way.

#
# NOTE: interesting conclusion, only Elastic Net - the only linear model extract from the model training process top ASVs specifically and no inflammatory marker.
# While, both Random Forest and XGBoost predicted mostly inflammatory marker and both models closely predicted on the same markers.
# I would conclude that it is better to compare RF and XGB as separate model, but keep Elastic Net alone, since it did not predict on any inflammatory markers.
#


#----------------------------------------------
# Receiver Operating Characteristic Comparison
#----------------------------------------------

# As mentioned before for EN, vector form is likely required to prevent errors;
ROC_EN <- roc(y_test_num, as.vector(EN_prediction))

# Plot all the AUC data for each model on the canvas
p35 <- plot(ROC_EN, col = "purple", main = "Comparison of ROC curves for each model")
p35 <- lines(ROC_RF, col = "red")
p35 <- lines(ROC_XGB, col = "blue")
p35 <- legend("bottomright", legend = c(paste("EN - AUC:", round(AUC_EN, 2)), paste("RF - AUC:", round(AUC_RF, 2)), paste("XGB - AUC:", round(AUC_XGB, 2))), col = c("purple", "red", "blue"), lwd = 2)

ggsave("plots/p34.png")
# Notice there are many immediate spikes from EN and overall XGBoost performed the best

# Although XGBoost performed the best it still reflect moderate predcitive power for this task.
# We can deduce that XGBoost has a 64% of predicting correctly whether a sample belongs to FV or UD group based on the complete ASV and inflammatory dataset.
# An AUC of this "low" score indicates that it is challenging to classify if a dietery treatment/intervention creates micrbial community shifts.

# The spikes in Elastic Net can be explained due to the linear nature of the model - the model is outputting a clustered probability threshold across smaller samples sizes.
# And Random forest is just guessing at this point, very unreliable


# Microbial shifts are extremely unique to each individual, therefore it could be that machine learning model return generally lower AUC because it cannot learn from direct patterns


#=====================================================
# Overview - what can be changed and further explored
#=====================================================

#------
# ASCA
#------
# 1) Retrieve the significant microbiome species of ASVs to see which are the key drivers -> extract plots and data frame > retrieve the loading

# 2) check if these said key drivers have any predictive power

#------------------
# Machine Learning
#------------------
# 3) check that the W0 ASVs and inflammatory markers do not change, else the data for the ones that change is insignificant and definitely carry no predictive power

# 4) Additionally, check what happens at W8 for ASVs to see the key changes


#==========================
# 1) RM-ASCA - key species
#==========================

# Retrieve all the information from the RM-ASCA pipeline from before into a dataframe
ASCA_loadings <- data.frame(Taxon = colnames(seq_asv_f_cleaned), Loading = ASCA_int_pca$rotation[, 1]) %>%
  # using the data that is filtered - f, outlier free - cleaned, and rCLR transformed
  mutate(Abs_Loading = abs(Loading)) %>% # the absolute value if needed so no negative for plotting, but it is interesting to see both 
  arrange(desc(Abs_Loading)) # arrange the data in descending order

# We shall plot the top 10 key species of taxa drivers;
p36 <- head(ASCA_loadings, 10) %>%
  ggplot(aes(x = reorder(Taxon, Loading), y = Loading, fill = Loading > 0)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values = c("TRUE" = "green", "FALSE" = "red")) + labs(title = "Top 10 Key taxonomic ASV drivers from RM-ASCA", subtitle = "Positive Loadings increase under FV, negatives decrease", x = "ASV / Taxon", y = "RM-ASCA Interaction Component 1 Loading") + theme(plot.title = element_text(hjust = 1), plot.subtitle = element_text(hjust = 1)) # this is to adjust the title location
ggsave("plots/p35.png")
# This displays a really nice plot of 6 positively loaded ASVs and 4 negatively.
# What this tells us - that ASV9 and ASV8 for example contribute with the most positive strength to the interaction effect (Treatment * time)

# Prevotellaceae prevotella - strict anaerobic and thrives in people consuming plant-rich diets that supply complex carbohydrates - could be related to FV
# Ruminococcaceae Faecalibacterium - anaerobic, digests fibers, they produce SCFAs (butyrate) which can regulate inflammation. It is said that this species decreases in elderly people and can be restored/boosted by eating pre-biotic diets or fermented dietary fibers.

# ASV36 and ASV28 contribute the most negative strength, the abundance of said ASVs are decreasing, not as driven in the study;

# Rikenellaceae Alistipes - Alistipes genus (when combined with the GI) provide protective effects against intestinal inflammation > specie specialised to the GI, 
# Ruminococcus torques group uncultured organism - bacterium detected in the gut microbiome through sequencing, and is not yet isolated, uncultured because it is difficult to grow certain species as they are sensitive to oxygen. they produce specific enzyme that breaks down mucin, the layer that protects your gut lining
#   by breaking mucin down, they produce further species that contribute overall to the micrbial health


# Now we can additionally, check the predictive power of said key ASVs from above - this can be done by training a random forest (works since there is no missing data)

#====================================
# 2) Predictive Power of Key Drivers
#====================================

if(!require("randomForest")) install.packages("randomForest")
library(randomForest)
packageVersion("randomForest") # Version: 4.7.1.2

# We looked at top 10 taxa, but not we can increase it a little to give the model more leeway to train on samples;
top_ASV_drivers <- head(ASCA_loadings, 50)$Taxon

target_binary <- ifelse(seq_meta_cleaned$Treatment == "FV", 1, 0) # gives binary outcomes to our target variable
unique_subjects <- unique(seq_meta_cleaned$ID) # same as before the extract the unique IDs, done again just for consistency for this section
# Define a train test split 80-20;
train_ids <- sample(unique_subjects, size = 0.8 * length(unique_subjects))

# create the train test split data for the model;
train_idx <- which(seq_meta_cleaned$ID %in% train_ids)
test_idx <- which(!seq_meta_cleaned$ID %in% train_ids)

asv_mat <- as.matrix(seq_asv_f_cleaned)

# Define the training and testing data based on the filtered, cleaned, and transformed ASVs
X_train <- asv_mat[train_idx, ]
X_test <- asv_mat[test_idx, ]

Y_train <- as.factor(ifelse(target_binary[train_idx] == 1, "FV", "UD")) # random forest cannot perform regression on numerical variables, so we need to factor them back, I tested it and got an error naturally
Y_test <- as.factor(ifelse(target_binary[test_idx] == 1, "FV", "UD")) # setting 1 to FV and 0 to UD

X_train_drivers <- as.data.frame(X_train)[, top_ASV_drivers] # now define the training and testing data sub-setted with the key taxonomic drivers from top 50
X_test_drivers <- as.data.frame(X_test)[, top_ASV_drivers]

# Run random forest on the selected top taxonomic species;
set.seed(1)
RF_drivers_model <- randomForest(x = X_train_drivers, y = Y_train, ntree = 500, importance = TRUE) # default values for now
RF_drivers_prediction <- predict(RF_drivers_model, X_test_drivers, type = "prob")[, "FV"] # validate on the test set on the target set for the intervention group

# Now we need the models metrics for performance and the predictive power for the ASVs;
ROC_drivers <- roc(ifelse(Y_test == "FV", 1, 0), RF_drivers_prediction) # extract the ROC of the model with 1 for FV and 0 for UD
print(paste("Key Taxa drivers only AUC:", round(auc(ROC_drivers), 2))) # round the ROC-AUC to 2 d.p. and print the result 
# 63% efficiency in predicting the Key taxonomic drivers to the intervention group FV, through randomForest

# To my surprise even with significantly driven taxa, the predictive performance on the model did not increase much -> AUC: 0.59 - 0.63;
# My thought process is that maybe the model is struggling to improve the predictive performance because it is taking in also the negatively affected taxa as key drivers?

#==================
# 3) Baseline ASVs
#==================

# Firstly, one must define the time points accordingly, from baseline (w0) to end point (w8);
w0_indices <- which(seq_meta_cleaned$time == "W0") # collects all W0 entries
w0_mat <- seq_asv_f_cleaned[w0_indices, ] # create a matrix sub setting the ASVs with only the W0 entries
w0_meta <- seq_meta_cleaned[w0_indices, ] # similarly creating a data item for the W0 meta data

# Starting with Bray-Curtis, since there are no negative values it should work and for Jaccard uses binary entries
# Again, we can use Bray-Curtis and Jaccard distance measures for the matrix to evaluate the baseline samples;

# Since both methods are faulty against negative values, Aitchison distance and Euclidean distances may be the best options
w0_distance_bray <- vegdist(w0_mat, method = "bray") # calculates the dissimilarity indices in the microbial community
w0_distance_jaccard <- vegdist(w0_mat, method = "jaccard")
# w0_distance_aitchison <- vegdist(w0_mat, method = "aitchison") !!! this throws an error as the method cannot be completed on CLR data !!!
w0_distance_euclidean <- vegdist(w0_mat, method = "euclidean") # seems euclidean distance is the best bet here

# https://github.com/vegandevs/vegan/issues/433, this discussion is about aitchison and euclidean distance measures pairing with CLR and rCLR, so we are on the right track.

# Again, running PERMANOVA to test the possibility of variation between samples over a set amount of iterations;
sample_stability <- adonis2(w0_distance_euclidean ~ ID, data = w0_meta, permutations = 999) # permutations = 999 is default
# The formula is the euclidean distance being calculated between each sample to see where the biggest variation lies
# this ran in a few seconds, while there are options for parallel computing so set the cores to compute the operation
print(sample_stability)

# Result - based on R^2 = ~0.88 (88%) of explained variance lies in this test. Therefore, we can deduce that a large proportion of variance is explained by the predictor variable, when looking at baseline samples
# Also, the p-value is significant - Pr(>F) = 0.001, p-value shows that the grouping pattern from the formula is significant that the sample IDs are clustered together at baseline, which is what we are looking for.

# This result tells me that the microbiome composition of subject 1 is 88% similar like subject 1 and the composition of subject 2 88% similar to subject 2, and so on.
# I can take this as significant that the microbiome composition of each individual is preserved to a high degree, but not completely 100%.

# Now since we have subject microbiome composition similarity, we can also check for specific subject composition differences at the end period of treatment.
# Investigate which species are differing from UD and FV at the end of the treatment and extract biological interpretation. It can signify which "rare" taxa the UD group are potentially missing, since they have not undergone the fermented intervention;


#=========================================
# 4) Differential Abundance Analysis - W8
#=========================================

# We can follow a similar pipeline to Differential Gene Expression Analysis (DGEA), but for microbiome - Differential Abundance Analysis (DAA)
# Here we can pin point which ASV species actually changed at the end-period of the treatment (Week 8), to distinguish the main differences between UD and FV.

# After further research, DAA method utilises ANOVA can also induce further bias and violates assumptions of this statistical testing because the samples are not independent!
# Therefore, a more robust method would be to use Linear Mixed Model (LMM) with an interaction term to see which ASVs changed with statistical significance over the course of the intervention period. Not just at the end of W8 since any ASV compositional "strength" that is larger than others, can inflate our results.

# The logic is that we want to test which ASVs perturbed over the "8 weeks" significantly enough between the two "groups", while keeping the subject "samples" as a random effect to account for their dependent nature.
# Notice that the words that are in quotations are the the variables that are encapsulated in our new interaction equation;

# https://cran.r-project.org/web/packages/lme4/index.html
# We will use lme4 package to complete LMM;
if(!require("lme4")) install.packages("lme4")
library(lme4)
packageVersion("lme4") # Version: 2.0.6

# This package allows us to test values from lme4 and automatically calculates p-values;
if(!require("lmerTest")) install.packages("lmerTest")
library(lmerTest)
packageVersion("lmerTest") # Version: 3.2.1

# initiate a results data frame to loop back all outputs;
perturbed_results <- data.frame()

# loop through all ASV columns;
ASV_species <- colnames(seq_asv_f_cleaned)[1:474] # there are 474 columns in the data frame we want to iterate through

# Additionally, for this Linear model analysis, the ASV table and meta data table need to be bind together, else the function below will give model = NULL as it will not be able to find the categorical variables;

LMM_combined <- seq_asv_f_cleaned %>%
  as.data.frame() %>%
  mutate(Treatment = seq_meta_cleaned$Treatment, time = seq_meta_cleaned$time, ID = seq_meta_cleaned$ID)

# sanity check for new added variables to be factors as usual;
class(LMM_combined$Treatment) # expect factors only 
class(LMM_combined$time)
class(LMM_combined$ID)

# create a for loop that a model will iterate through each ASV and test the perturbed significance;
for(ASV in ASV_species) {
  # fit the Linear Mixed Model
  LMM_formula <- as.formula(paste(ASV, "~ Treatment * time + (1 | ID)")) # the interaction formula we mentioned before
  LMM_model <- tryCatch(lmer(LMM_formula, data = LMM_combined), error = function(e) NULL) # tryCatch runs the first entry of the bracket and if there is an error - throws NULL just so that running of the loop doesn't stop abruptly
  
  if(!is.null(LMM_model)) { # if the Linear Mixed Model is NOT null then;
    LMM_model_summary <- summary(LMM_model)$coefficients # create a summary of the model statistical metrics - need the coefficients and p-values to gauge the significance
    
    # Since our focus here is the interaction formula, specifically; Treatment * time;
    if("TreatmentFV:timeW8" %in% rownames(LMM_model_summary)) { # from the summary table both intervention group and time Week 8 is extracted
      interaction_p_value <- LMM_model_summary["TreatmentFV:timeW8", "Pr(>|t|)"] # collect the significant intervention group summary statistics at Week 8
      LMM_estimate_summary <- LMM_model_summary["TreatmentFV:timeW8", "Estimate"] # gather the mode estimates at the end of the period for intervention
      perturbed_results <- rbind(perturbed_results, data.frame(ASV = ASV, Estimate = LMM_estimate_summary, P_value = interaction_p_value)) # bind all results to a single data object
    }
  }
}
print(head(perturbed_results))

# NOTE: be careful with the variable names created automatically with the lme4 packages i.e. LMM_model and LMM_model_summary, usually NULL is output when the variable names used in quotations are incorrectly spelt - syntax

# Moreover, since this for loop is an iteration over all 474 ASVs one by one, there is a high chance for false positives from statistical chance.
# To account for this one can perform Benjamini-Hochberg correction method for False Detection Rate (FDR)

# Utilising the p-adjust function to adjust the p-values accordingly;
# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust

perturbed_results$Adjusted_P_value <- p.adjust(perturbed_results$P_value, method = "BH") 
# And we set a p-value significance threshold of 0.05 to filter out unnecessary data;
significant_perturbed_results <- subset(perturbed_results, Adjusted_P_value < 0.05)
print(significant_perturbed_results)

# There are no significant ASVs that changed between the two groups of intervention and control over the 8 Weeks ACCORDING to the threshold of an adjusted p-value below 0.05

# We can still experimentally check without the Benjamini-Hochberg correction for FDR, but then the results cannot be taken in complete faith since the result now can signify potential false positives (PFP);
perturbed_results_PFP <- subset(perturbed_results, P_value < 0.05)
print(perturbed_results_PFP)

nrow(perturbed_results_PFP) / nrow(perturbed_results) * 100 # the percentage of all perturbed Linear Mixed Model ASVs that are significant, but can carry a false positive pattern.
# Just under 7% (6.96% = 33 ASVs) would be classed as significant of all 474 cleaned ASVs

perturbed_results_PFP_top_10 <- perturbed_results_PFP %>% # get the significant ASVs
  mutate(Abs_Estimate = abs(Estimate)) %>% # use the beta diversity estimate coefficient to gather the ASVs that were perturbed the most over the 8 weeks
  slice_head(n = 10) %>% # gather top 10
  mutate(Direction = ifelse(Estimate > 0, "Increased (FV)", "Decreased (FV)")) # set a direction that can be seen on the plot: Positive - so the ASV perturbation increased within the intervention group and Negative - ASV perturbation decreased within the intervention group

perturbed_results_PFP_top_10$ASV <- reorder(perturbed_results_PFP_top_10$ASV, perturbed_results_PFP_top_10$Estimate) # reorders the values so that the largest at placed at the top for the Estimates and ASVs - nicer for the plot visualisation

p37 <- ggplot(perturbed_results_PFP_top_10, aes(x = ASV, y = Estimate, fill = Direction)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values = c("Increased (FV)" = "green", "Decreased (FV)" = "red")) + labs(title = "Top 10 ASVs of the Intervention Group (W8) - PFP (P value < 0.5)", x = "ASVs" , y = "Estimate (Beta diversity)") + theme(plot.title = element_text(hjust = 0.8))
ggsave("plots/p36.png")
# Interestingly enough, the results are nearly the opposite from RM-ASCA pipeline in terms of positive and negative contribution of species.
# Additionally, while again stressing the fact that there is a possibility of false positives, but this plot showcases the ASVs that contributed the most strength - perturbed the most over the 8 weeks of intervention compared to the control group
# So eat more fermented vegetables, so you can increase your microbial diversity and be beta diversity superior to others.
