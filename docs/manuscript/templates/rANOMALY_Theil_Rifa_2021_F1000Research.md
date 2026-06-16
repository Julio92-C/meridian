# rANOMALY: AmplicoN wOrkflow for Microbial community AnaLYsis

**PMID:** 33537122  
**DOI:** [10.12688/f1000research.27268.1](https://doi.org/10.12688/f1000research.27268.1)  
**PMC:** PMC7836088  

Source: PubMed Central full-text retrieval, 2026-05-29.

---

Introduction

Studies of microbial communities tends to become a daily routine analysis for lots of laboratories and the main method to explore microbial diversity is metabarcoding, which is an amplicon targeted sequencing method (16S for bacteria and ITS for fungi). Metabarcoding generates a large amount of data and a lot of applications already exist for their processing (FROGS, qiime). Methods and software are continuously evolving and the main challenge for bioinformaticians is to implement the most recent and effective ones in their analysis. Here we present rANOMALY, a scalable and lightweight R package which is able to handle every step of a metabarcoding analysis, from read cleaning, contaminant filtering, taxonomic assignment, to advanced statistical analysis. rANOMALY is fully implemented in R language in which each step correspond to one function, allowing to easy implementation of new features or tools while being easy to use and maintain. The package allows the workflow to be executed on any R environment. rANOMALY only needs a CSV table describing the metadata for each sample, and a folder containing the corresponding fastq files as input. It can produce high quality figures for Rmarkdown reports along with statistical tests ready for publication.

The workflow is illustrated in.

Figure 1.

rANOMALY workflow.

Implementation

rANOMALY is an R package depending on other CRAN, Bioconductor, and git R packages. It is easy to install via devtools::install_git function.

Raw sequence processing

Samples must have been previously demultiplexed into one file per sample with the file name following this syntax:. The denoising process is handled using theR packagewhich produces amplicon sequence variants (ASV) as a taxonomic unit. This improves resolution of the potential presence of microbial organisms by using a prediction model to correct sequencing errors before aggregating similar sequences. rANOMALY handles processing any region of 16S (V1 to V9) and ITS amplicon (ITS1, ITS2) sequences, in which an additional step ofremoves ITS probes left in some short sequences. For 16S amplicon, primers are trimmed based on the primer sequence length. ASV identifiers are sequences translated into MD5 hashes which are unique identifiers based on the DNA sequences offering the possibility to be compared between projects. This step results in an object with representing sequences, a raw ASV counts table and a text file containing statistics from the denoising process.

Taxonomic assignment

Taxonomic assignments of ASVs are carried out by IDTAXA (part of thepackage), an algorithm based on a machine learning method. We implemented a functionality to compute assignment with two reference databases. For example, 16S amplicon sequences can be assigned with an environment specific database (DAIRYdb, HITdb, MIDAS) and a general database like SILVAor GreenGenes. The assignment with the best confidence and the lowest rank is kept, thus increasing assignment depth and accuracy. As taxonomy can differ between reference databases, we implemented a taxonomy validation step in rANOMALY to unravel taxonomy inconsistency like taxa with multiple ancestors or empty ranks. As supplementary features, rANOMALY functions allow users to create their own IDTAXA formatted specific reference database. The first step consists of filling the taxonomic table empty fields with the last known rank, and checking for taxonomy incongruencies as in the assignment function. A taxid file as used with RDP classifieris constructed and then a last function takes as input, the corrected taxonomy table, the fasta file and the taxid file to generate the IDTAXA formatted database.

Phylogenetic tree

Phylogenetic tree is generated in three steps. First, sequences are aligned with AlignSeqs function from DECIPHER packageby the guided-tree method. Then, distance matrix is calculated with the dist.ml function from thepackage. And finally, neighbour joining and pml function computes the likelihood of the phylogenetic tree.

The abundance, metadata, taxonomic table, reference sequences and phylogenetic tree are merged into aobject.

Decontamination

ASVs have the advantage of enabling the distinction between contaminants and the real community. We have integrated thepackageinto the workflow. Indeed using OTU based clustering methods can agglomerate contaminant sequences with real sample sequences, the whole cluster could hence be considered as contaminant by mistake. Working with ASVs allows the use of R package decontam which will sensitively exclude contaminant ASVs. It integrates two main methods, one based on the prevalence of the contaminant in the control samples, and another one based on the DNA concentration of the samples. Moreover, our decontamination step allows users to apply various filters such as low ASVs frequency, low ASVs prevalence in real samples, and the minimum number of reads per sample.

Graphical plots and statistical analyses

Statistical analyses are key features of rANOMALY workflow. Main descriptive analyses are integrated thanks to phyloseq functionalities. In addition, we automatized graphical representations and advanced statistical tests for alpha, beta diversity and composition plots. Above all, we have included the four most up-to-date differential analyses to assess differentially abundant features.

Community composition plots, alpha and beta diversity analyses

rANOMALY allows users to explore microbial community composition with three different types of plots : classical interactive bar plotof raw and relative taxa abundances, rarefaction curves to check sampling effort, and Krona interactive pie charts.

indices representing the specific richness are calculated (Richness, Simpson, Shannon...) with theR package. We added statistical tests such as multi-factors analysis of variance and pairwise Wilcoxon tests to assess significant differences between tested categories. We included repeated measures ANOVA to handle within-subjects variation. For example, it can be used when there are measures taken on the same individuals at different time points. This step outputs graphical representations and tables with results of statistical tests which are both saved in files and returned as list objects for markdown reports.

analysis allows users to estimate the community differences between two samples, it is also based on the vegan package. rANOMALY can calculate all different distances such as the one based on ASVs abundance (BrayCurtis), rank based Jaccard indices, and phylogenetic distances as UniFrac, weighted Unifrac. Graphical representations PCoA, NMDS and more are available. This analysis can be processed at different taxonomic levels and categorical factors chosen by the user. The additional statistical test of PERMANOVA uses the distance matrix to determine if microbial communities of sample groups are significantly different from each other. We added a pairwise PERMANOVAtest to confirm significant differences between specific group of samples.

Differential analysis

Differential analyses are meant to assess potential differentially abundant taxas between tested conditions chosen by users. rANOMALY wraps three methods:

The use of multiple methods for differential analysis allows the user to investigate which feature can be considered as differentially abundant between conditions. rANOMALY function using DESeq2 and metagenomeSeq outputs tables and plots with significant features. Metacoder outputs heat tree plots allowing users to infer differentially abundant features at each taxonomy rank and their position on the phylogenetic tree. With basic data management, and according to the identification of all significant differentially features (i.e. ASVs, genus, species) from the three methods, we generate a single table to find significant features in one or more methods and ease the interpretation. Additional information are added to the final table, like mean relative abundance for each feature and condition. A column in which condition features are significantly more abundant, features taxonomies and related sequences are added to the aggregated table. To complete differential analysis, we have included the well recognized PLS-DA from thepackage. It is a supervised classification method allowing users to identify features discriminating the sample groups.

uses negative binomial generalized linear model with a variance stabilizing transformation on abundances.

uses a zero-inflated log-normal model with cumulative sum scaling normalization.

applies a total sample sum normalisation and uses a non-parametric Wilcoxon Rank Sum test to compare the log ratio of mean proportions.

Additional features and export

Functions and procedures are available to help the user to generate additional figures or to export the data to third party software. For instance, shared taxa between conditions are useful to explore. We use a function that can generate Venn diagrams, or for more complex visualisation, we can produce files readable by Cytoscapeto produce shared taxa networks. Krona diagrams can be displayed to explore sample microbial composition, where samples can be merged by a specified factor. rANOMALY allows users to generate inputs for, which is a graphical software that provides statistical hypothesis tests and exploratory plots for analysing taxonomic profiles.

Operation

rANOMALY requires R 3.6.3 or upper and can be run on any operating system with common specifications (1Go disk space, 4Go RAM, multicore CPU is recommended).

Use case

For this example we are using a dataset from Fretinstudyin which samples are from four different environments: cow milk, cow cheese (rind and core) and cow teat skin. This dataset and metadata are available on NCBI-SRA website: BioProject accession. To ease access to this dataset, fastq files along with pre-formatted metadata are available on this.

Install

Up to date code is hosted by the, users can simply download and install the package in R console with following command lines.

()()()()

Processing of raw sequences

The first step is to define ASVs thanks to thepackage. In, only one function is needed to compute all the different steps require from this package. Here sample names are extracted from the file name, thus be sure that files name match samples name in the metadata file.

(=,=,=,=,=,=,=,=,=)

Main outputs of this function are:

()(())

function uses IDTAXA function from DECIPHER package, and allows to use two different databases. It keeps the best assignment on two criteria, resolution (depth in taxonomy assignment) and confidence (value givenby IDTAXA). The final taxonomy is validated by checking for multiple ancestors (same species assigned to different genus) and incongruities (empty fields or incomplete lineage) correction step.

We share the latest databases we use in the IDTAXA format in this. Users can also generate your own IDTAXA
formatted database following those instructions and scripts we provide at.

(=(,),=,=)

Main file outputs:

The phylogenetic tree from the representative sequences is generated usingandpackages.

tree

= generate_tree_fun

(dada_res)

To create aobject, we need to merge four objects and one file:

(,)

Thefunction usesR package along with control samples (PCR control) to filter out contaminants. Thepackage offers two main methods, frequency and prevalence (users can also combine those methods). For frequency method, it is mandatory to have the DNA concentration of each sample in the phyloseq object (and hence in the). The prevalence method does not need DNA quantification, this method allows to compare presence/absence of ASV between real samples and control samples and then identify contaminants.

sequencing plateforms often quantify the DNA before sequencing, but do not usually give the information. Just ask for it ;).

Our function integrates the basic ASV frequencyand minimum prevalence in overall samples filtering. We have also included an option to filter out ASV based on their taxa names for known laboratory recurrent contaminants.

Main outputs:

Here we are going to filter out ASVs representing less than 0.1% () of the reads and that are present in more than 4 samples (). Moreover, we are excluding "unassigned" taxas for this use case. Our sample dataset do not contain control samples, this step will be skipped.

(=,=,=,,,,,)

We obtain the final phyloseq object used for downstream analysis:

> data_filtered
phyloseq-class experiment-level object

otu_table

()

OTU Table:         [ 59 taxa and 48 samples ]

sample_data

()

Sample Data:       [ 48 samples by 34 sample variables ]

tax_table

()

Taxonomy Table:    [ 59 taxa by 7 taxonomic ranks ]

phy_tree

()

Phylogenetic Tree: [ 59 tips and 57 internal nodes ]

refseq

()

DNAStringSet:      [ 59 reference sequences ]

summarizes the read number after each filtering step ().

the raw ASV table.

fasta file with all representative sequences for each ASV.

with savedlist containing raw ASV table and representative sequences in objects.

with taxonomy in phyloseq format inobject.

the final assignation table () outputed in CSV format.

raw assignations from the two databases, mainly for debugging.

the ASV tableand the representative sequencesinvariable.

a taxonomy table ().

the phylogenetic tree ().

metadata from from csv file.

with contaminant filtered phyloseq object named.

list of filtered ASVs for each filtering step.

Krona plot before and after filtering.

and

venn diagram showing the repartition of filtered ASVs by decontamination methods.

Table 1.

Read tracking table overview.

sample.id

input

ï¬ltered

denoisedF

nonchim

SRR6365127.fastq

35690

33511

33402

33402

SRR6365128.fastq

38871

38183

38111

38046

SRR6365129.fastq

50970

49765

49631

49092

SRR6365130.fastq

44207

42657

42524

42524

SRR6365131.fastq

68008

66517

66257

66087

SRR6365132.fastq

44014

42158

41945

41848

SRR6365133.fastq

36704

35884

35799

35444

...

...

...

...

...

Plots, diversity and statistic analyses

In order to observe the sampling depth of each sample we start by plotting rarefaction curves. Those plots are generated bywhich makes them interactive. ()

,)

Thefunction allows user to generate interactive community composition plot.presents the composition plot with relative abundances for the top 20 genera existing in our samples. The function allows to plot at different taxonomy rank and to modify the number of taxa to show.

Here two arguments controlling the composition plot aesthetics:

factor shows very different bacterial community between milk, cheese and cow teats environments.

,=,=,=,=,=)

Thefunction can computes various alpha diversity indexes. It uses thefunction from(Available measures : Observed, Chao1, ACE, Shannon, Simpson, InvSimpson, Fisher). Here we calculate ASV richness and Shannon index and carry out an analysis of variance on thefactor. Sequencing depth is automatically taken into account in this test. A pairwise wilcoxon test is added to ANOVA to define which group might be significantly different from others.shows boxplots of diversity indices, cow teats environment has much more ASV than other environments. Shannon index reveals more differences between cheese rinds and cheese cores. Cheese rinds show a higher Shannon index highlighting a more balanced bacterial community.

(=,,(,))

Results of the analysis of variance and pairwise wilcoxon test on Shannon index:

# > divAlpha$Shannon
# $anova
#                  Df Sum Sq Mean Sq F value  Pr(>F)
# Depth             1   0.83    0.83   9.203 0.00409 **
# source_location  3  98.37   32.79 365.661 < 2e-16 ***
# Residuals        43   3.86    0.09
# ---
# Signif. codes:  0 â€˜***â€™ 0.001 â€˜**â€™ 0.01 â€˜*â€™ 0.05 â€˜.â€™ 0.1 â€˜ â€™ 1
#
# $wilcox_col1
#           milk past_rind raw_rind
# past_rind    0        NA       NA
# raw_rind     0         0       NA
# teat         0         0        0

Thefunction returns a list which contains:

And for each of the computed indices :

Thefunction allows to generate specific tests and figures ready to publish in rmarkdown report as in the example below. It is based on the vegan package functionfor the distance calculation andin addition tofuntion for the ordination plot.

We include statistical tests to ease the interpretation of results. A permutational ANOVA is carried out on matrix distance to compare groups by testing if centroids and dispersion are equivalent for all groups. User have to informargument and optionally(covariable) to assess PERMANOVA to determine significant differences between groups. A pairwise-PERMANOVA is processed to determine which condition is significantly different from another (based on p-value).

As a return, you will get a list that contains:

Here we present results of beta diversity analysis onfactor,show ordination plot and it confirms the big differences between community of cheese, milk and cow teats.

(,,,,

The permanova tests on BrayCurtis distance shows significant p-value for thefactor.is used to define which level of the factor tested is significantly different from others.

divBeta

$

permanova

# Permutation: free
# Number of permutations: 1000
#
# Terms added sequentially (first to last)
#
#                 Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)
# Depth            1    2.4044  2.4044  27.072 0.14822 0.000999 ***
# source_location  3    9.9981  3.3327  37.523 0.61634 0.000999 ***
# Residuals       43    3.8191  0.0888         0.23544
# Total           47   16.2217                 1.00000
# ---
# Signif. codes:  0 â€˜***â€™ 0.001 â€˜**â€™ 0.01 â€˜*â€™ 0.05 â€˜.â€™ 0.1 â€˜ â€™ 1

divBeta

$

pairwisepermanova

#                        pairs Df SumsOfSqs   F.Model        R2 p.value p.adjusted sig
# 1        milk vs cheese_core  1 5.2215038 85.356493 0.7950753   0.001     0.0012   *
# 2        milk vs cheese_rind  1 4.8836384 53.310404 0.7078757   0.001     0.0012   *
# 3               milk vs teat  1 3.9377389 26.676358 0.5480352   0.001     0.0012   *
# 4 cheese_core vs cheese_rind  1 0.2640387  4.573476 0.1721068   0.008     0.0080   *
# 5        cheese_core vs teat  1 4.7206284 41.504931 0.6535702   0.001     0.0012   *
# 6        cheese_rind vs teat  1 4.3806123 30.384779 0.5800307   0.001     0.0012   *

We choose three different methods to process differential analysis which is a key step of the workflow. The main advantage of the use of multiple methods is to cross validate deferentially abundant taxa between tested conditions. For this use case, we choose to focus on milk and cow teat environment to compare community at genus level.

Metacoder is the most simple differential analysis tool of the three. Counts are normalized by total sum scaling to minimize the sample sequencing depth effect and it uses a Kruskal-Wallis test to determine significant differences between groups. Thefunction allows the user to choose the taxonomic, which factor to the test (), and a specific pairwise comparison () to launch the differential analysis.

It produces pretty graphical trees, representing taxas present in both groups and coloring branches depending on which group this taxa is more abundant (). Two trees are produced, a raw one, where everything is displayed and a filtered one where only significant features are represented (p-value <= 0.05).

It produces pretty graphical trees, representing taxas present in both groups and coloring branches depending on which group this taxa is more abundant (). Two trees are produced, a raw one, where everything is displayed and a filtered one where only significant features are represented (p-value <= 0.05).

(,=,=,=,=,=,=)

Main output is a list with :

DESeq2 is a widely used method, primarily for RNAseq applications, for assessing differentially expressed genes between controlled conditions. Its use for metabarcoding datas is sensibly the same and well documented. Theallows to process differential analysis as, and users can choose the taxonomic rank, the factor to test and which condition to compare. DESeq2 algorithm uses negative binomial generalized linear models with VST normalization (Variance Stabilizing Transformation).

Main output is a list with:

(=,,,,)

MetagenomeSeq uses a normalization method able to control for biases in measurements across taxonomic features and a mixture model that implements a zero-inflated Gaussian distribution to account for varying depths of coverage. Asreturns a table with statistics and a plot with significant features for each comparison.

(=,,,,)

Thefunction allows to merge the results from the three differential analyses methods computed previously to obtain one unique table with all informations of significant differentially abundant features.shows the most significant and differential abundant Genera between the two environments.

,,,,,,,,)

The generated table include the following fields:

Here is an overview of thetable informing in which methods each feature is significant, their DESeq2 LogFoldChange value, taxonomy and representative sequences:

,,

Figure 2.

Rarefaction plot for all samples, separated by source_location factor.

Figure 3.

Composition plot of relative abundance of 20 most abundant genus.

Figure 4.

Boxplot of diversity alpha indexes.

Figure 5.

Ordination plot: MDS on Braycurtis distance.

Figure 6.

Heattree of significant features generated by metacoder package.

Figure 7.

Top significant and differentially abundant genera between Milk and Teat samples.

option order the sample along the X axis.

option control labels of the X axis.if user donâ€™t want the sample to be renamed.

boxplots comparing conditions with chosen indices. ()

a table of indices values. ()

an ANOVA analysis. ()

a pairwise wilcox test result comparing conditions and giving the pvalue of each comparison tested:wilcox test results on the first or unique factor,: wilcox test results on the second factor,: wilcox test results on collapsed factor 1 and factor 2.

a mixture model if your dataset includes repeated measures, ie.option. ()

An ordination plot ().

The permANOVA results ().

The pairwise permANOVA ()

for each comparison (), two heattree, one with all features () and an other one with only significant features ().

a table with all wilcoxon test results ().

a plot showing Log2FoldChange value of each significant feature ().

a table with statistics (LogFoldChange, pvalue, adjusted pvalue...) ().

seqid: ASV ID.

Comparaison: Tested comparison.

Deseq / metagenomeSeq / metacoder: differentially abundant with this method (0 no or 1 yes).

sumMethods: sum of methods in which feature is significant.

DESeqLFC: Log Fold Change value as calculated in DESeq2.

absDESeqLFC: absolute value of Log Fold Change value as calculated in DESeq2.

MeanRelAbcond1 / MeanRelAbcond2: Means relative abundance in condition 1 and 2.

Condition: in which the mean feature relative abundance is higher.

Taxonomy and representative sequence.

Miscellaneous function

User can generate an interactive heatmap withfunction to explore relative abundance of top taxa through samples as showed in.

(=,=,,)

shows Venn diagram comparing the three environment of this study, this method is useful to determine shared taxa between group of samples. This function usesfunction which handles up to 7 groups comparison.

,,)

generate a table allowing user to find which taxa are shared between conditions.

,,)))

Figure 8.

Heatmap of top genuses relative abundance.

Figure 9.

Venn plot of shared genuses.

Export and compatibility

allows to generate input files (SIF format) forwhich is useful to visualise shared taxa and easily modify each nodes and arrows position and aesthetic.

function allows user to generate tabulated format abundance table at different taxonomic rank with following commnands:

(,,,(),}

function creates two files that can be imported into the STAMP software.

(=,=)

function allows user to import data from other bioinformatic pipeline like FROGS, Qiime2. This function needs tabulated ASV, taxonomy, metadata and DNA sequence table as inputs. It can generate phylogenetic tree if missing and output a phyloseq object ready for downstream analyses.

Conclusions

rANOMALY allows users to handle metagenomic data from raw sequences quality control to final differential analysis with ready to publish results in an easy and reproducible manner. Users have access to all sources of the rANOMALY package that can be deployed on any operating system or server allowing them to analyse anything from a few samples to several thousand. This workflow combines all of the latest developments in the field: the use of high resolution amplicon sequence variant, contaminant filtering, double automatic taxonomic assignation, integrated statistical analyses and four differential analyses with cross validation. rANOMALY help users those who donâ€™t have big programming skills. And since rANOMALY uses phyloseq objects and standards, original functions from the phyloseq package can still be used to split, filter and select specific samples previously to be visualized and/or tested with rANOMALY functions. For exploratory analysis and interactive experience, a shiny application based on this workflow is under active development and will be shared with the community.

Software availability

Up-to-date source code, and tutorials are available at:.

Package documentation is also provided at:

Archived source code as at time of publication are available from:

License:

Here we show all packages used in this workflow and their version numbers:

sessionInfo

()

# R version 3.6.3 (2020-02-29)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.5 LTS
#
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
#
# locale:
#  [1] LC_CTYPE=fr_FR.UTF-8       LC_NUMERIC=C               LC_TIME=fr_FR.UTF-8        LC_COLLATE=fr_FR.UTF-8
#  [5] LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=fr_FR.UTF-8    LC_PAPER=fr_FR.UTF-8       LC_NAME=C
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C
#
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
# [1] ranomaly_0.0.0.9000
#
# loaded via a namespace (and not attached):
#   [1] taxa_0.3.4                  tidyselect_1.1.0            htmlwidgets_1.5.1          RSQLite_2.2.0
#   [5] AnnotationDbi_1.48.0        grid_3.6.3                  BiocParallel_1.20.1        Rtsne_0.15
#   [9] devtools_2.3.2              IHW_1.14.0                  munsell_0.5.0              codetools_0.2-16
#  [13] withr_2.2.0                 colorspace_1.4-1            Biobase_2.46.0             phyloseq_1.30.0
#  [17] knitr_1.30                  rstudioapi_0.11             stats4_3.6.3               ggsignif_0.6.0
#  [21] slam_0.1-47                 GenomeInfoDbData_1.2.2      lpsymphony_1.14.0          hwriter_1.3.2
#  [25] bit64_0.9-7.1               rhdf5_2.30.1                rprojroot_1.3-2            vctrs_0.3.2
#  [29] generics_0.0.2              xfun_0.16                   lambda.r_1.2.4             R6_2.4.1
#  [33] GenomeInfoDb_1.22.1         locfit_1.5-9.4              bitops_1.0-6               microbiome_1.8.0
#  [37] DelayedArray_0.12.3         assertthat_0.2.1            scales_1.1.1               gtable_0.3.0
#  [41] processx_3.4.4              phangorn_2.5.5              rlang_0.4.7                genefilter_1.68.0
#  [45] splines_3.6.3               lazyeval_0.2.2              rstatix_0.6.0              broom_0.7.0
#  [49] reshape2_1.4.4              abind_1.4-5                 backports_1.1.8            tools_3.6.3
#  [53] usethis_1.6.3               ggplot2_3.3.2               ellipsis_0.3.1             gplots_3.1.0
#  [57] decontam_1.6.0              biomformat_1.14.0           RColorBrewer_1.1-2         BiocGenerics_0.32.0
#  [61] sessioninfo_1.1.1           Rcpp_1.0.5                  plyr_1.8.6                 zlibbioc_1.32.0
#  [65] psadd_0.1.3                 purrr_0.3.4                 RCurl_1.98-1.2             ps_1.4.0
#  [69] prettyunits_1.1.1           ggpubr_0.4.0                Wrench_1.4.0               viridis_0.5.1
#  [73] S4Vectors_0.24.4            SummarizedExperiment_1.16.1 haven_2.3.1                cluster_2.1.0
#  [77] fs_1.4.2                    DECIPHER_2.14.0             magrittr_1.5               RSpectra_0.16-0
#  [81] data.table_1.13.0           futile.options_1.0.1        pairwiseAdonis_0.0.1       openxlsx_4.1.5
#  [85] ranacapa_0.1.0              matrixStats_0.56.0          pkgload_1.1.0              evaluate_0.14
#  [89] hms_0.5.3                   xtable_1.8-4                XML_4.0-0                  VennDiagram_1.6.20
#  [93] rio_0.5.16                  jpeg_0.1-8.1                readxl_1.3.1               IRanges_2.20.2
#  [97] gridExtra_2.3               shape_1.4.4                 testthat_3.0.0             compiler_3.6.3
# [101] ellipse_0.4.2               tibble_3.0.3                KernSmooth_2.23-17         crayon_1.3.4
# [105] htmltools_0.5.0             venn_1.9                    corpcor_1.6.9              mgcv_1.8-33
# [109] tidyr_1.1.0                 geneplotter_1.64.0          RcppParallel_5.0.2         DBI_1.1.0
# [113] formatR_1.7                 MASS_7.3-53                 ShortRead_1.44.3           Matrix_1.2-18
# [117] ade4_1.7-15                 car_3.0-10                  permute_0.9-5              cli_2.0.2
# [121] quadprog_1.5-8              parallel_3.6.3              igraph_1.2.5               GenomicRanges_1.38.0
# [125] forcats_0.5.0               pkgconfig_2.0.3             GenomicAlignments_1.22.1   foreign_0.8-76
# [129] plotly_4.9.2.1              foreach_1.5.0               rARPACK_0.11-0             annotate_1.64.0
# [133] admisc_0.9                  multtest_2.42.0             XVector_0.26.0             stringr_1.4.0
# [137] callr_3.5.1                 digest_0.6.25               vegan_2.5-6                dada2_1.14.1
# [141] Biostrings_2.54.0           rmarkdown_2.5               cellranger_1.1.0           fastmatch_1.1-0
# [145] curl_4.3                    Rsamtools_2.2.3             gtools_3.8.2               lifecycle_0.2.0
# [149] nlme_3.1-149                jsonlite_1.7.1              Rhdf5lib_1.8.0             mixOmics_6.10.9
# [153] carData_3.0-4               futile.logger_1.4.3         desc_1.2.0                 viridisLite_0.3.0
# [157] limma_3.42.2                fansi_0.4.1                 pillar_1.4.6               metacoder_0.3.4
# [161] lattice_0.20-41             httr_1.4.2                  plotrix_3.7-8              pkgbuild_1.1.0
# [165] survival_3.1-12             glue_1.4.1                  remotes_2.2.0              zip_2.0.4
# [169] fdrtool_1.2.15              png_0.1-7                   iterators_1.0.12           glmnet_4.0-2
# [173] bit_1.1-15.2                stringi_1.4.6               metagenomeSeq_1.28.2       blob_1.2.1
# [177] DESeq2_1.29.13              latticeExtra_0.6-29         caTools_1.18.0             memoise_1.1.0
# [181] dplyr_1.0.0                 ape_5.4
