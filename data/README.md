# Effects of Fermented Vegetables on the Gut Microbiota for Prevention of Cardiovascular Disease

[Access this dataset on Dryad](https://doi.org/10.5061/dryad.547d7wmfk)
This dataset includes data on amplicon sequence variants (ASVs) derived from next-generation sequencing of the V4 region of 16S rRNA isolated from human stool samples. Sequencing was performed in a pair-end modality on the Illumina MiSeq 500 platform rendering 2 x 150 bp paired-end sequences. These data were obtained from subjects who were either randomized into a fermented vegetable group, consuming 100 g of fermented vegetables on at least 5 days per week for 8 weeks, or a control group who were asked to not change anything in their diet. All subjects had at least one risk factor for cardiovascular disease and ages ranged from 35 to 60 years of age.

## Description of the data and file structure

#### File name: ASV\_table\_dryad.xlsx

File description: This file contains the output from next-generation sequencing using the Illumina MiSeq platform. It contains data on microbial alpha and beta diversity, as well as amplicon sequencing variants (ASVs).

The file includes the following data:

* Column A: subject ID
* Column B: Timepoint – 1 is for baseline and 2 is for follow up (8 weeks later)
* Column C: Treatment – 1 = control and 2 = fermented vegetable
* Column D: Observed features – number of ASVs detected
* Column E: Chao1 – indicator of species richness (total number of species in a sample)
* Column F: Simpson – indicator of species evenness (proportional distribution of the number of each species in a sample)
* Columns G, H, I: Bray Curtis index PCOA plot distances
* Columns J, K, L: Jaccard index PCOA plot distances
* Columns M, N, O: Weighted unifrac PCOA plot distances
* Columns P, Q, R: Unweighted unifrac PCOA plot distances
* All other columns: Amplicon sequence variants (ASVs)

Subjects in this study collected three stool samples at baseline and three stool samples at 8 weeks, therefore, there are six rows per person.

#### File name: Inflammatory\_markers

File description: This file contains data on inflammatory markers measured in serum collected from study participants at baseline  (time point 1) and at week 8 (time point 2). 

The file includes the following data:

* Column A: subject ID
* Column B: TMAO_1ngmL - Serum Trimethylamine oxide (ng/mL) at time point 1
* Column C: TMAO_2ngmL - Serum Trimethylamine oxide (ng/mL) at time point 2 (8 weeks later)
* Column D: LBP1_ngmL - Serum Lipopolysaccharide binding protein (ng/mL) at time point 1
* Column E: LBP2_ngmL - Serum Lipopolysaccharide binding protein (ng/mL) at time point 2 (8 weeks later)
* Column F: Ferr1_ngmL - Serum Ferritin (ng/mL) at time point 1
* Column G: Ferr2_ngmL - Serum Ferritin (ng/mL) at time point 2 (8 weeks later)
* Column H: CRP1_ngmL - Serum C-reactive protein (ng/mL) at time point 1
* Column I: CRP2_ngmL - Serum C-reactive protein (ng/mL) at time point 2 (8 weeks later)
* Column J: LOX1_pgmL - Serum Oxidized LDL Receptor (pg/mL) at time point 1
* Column K: LOX2_pgmL - Serum Oxidized LDL Receptor (pg/mL) at time point 2 (8 weeks later)
* Column L: Fiaf1_ngmL - Serum angiopoietin-like protein 4 at time point 1
* Column M: Fiaf2_ngmL - Serum angiopoietin-like protein 4 at time point 2 (8 weeks later)

Cells marked with a '.' indicate a missing value 

