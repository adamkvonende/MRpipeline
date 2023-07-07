
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Hopewell Group MR pipeline

:information_source: This pipeline is still under active development.  
Check the [NEWS](NEWS.md) to learn more about what has been
updated/modified!

## Overview

This pipeline is designed to facilitate genome-wide and cis-Mendelian
randomisation (MR) analyses. In addition to applying various univariable
MR methods, the pipeline includes options for data pre-processing,
harmonisation, and clumping. The function generates results tables,
relevant statistics (e.g. heterogeneity statistics), and a forest plot
of the results, as well as additional tables and information to
facilitate further analyses.

## Installation

The pipeline relies on one major function:

- **`perform_mr()`**

The easiest way to get started is to use the [Example](#example)
provided below.

The functions relies on various R packages. Some of these
([`rio`](https://cran.r-project.org/web/packages/rio/index.html),
[`tidyverse`](https://cran.r-project.org/web/packages/tidyverse/index.html),
[`data.table`](https://cran.r-project.org/web/packages/data.table/index.html),
[`MendelianRandomization`](https://cran.r-project.org/web/packages/MendelianRandomization/index.html),
[`testthat`](https://cran.r-project.org/web/packages/testthat/index.html),
[`grid`](https://cran.r-project.org/web/packages/grid/index.html),
[`gridExtra`](https://cran.r-project.org/web/packages/gridExtra/index.html),
[`ieugwasr`](https://cran.r-project.org/web/packages/ieugwasr/index.html),
[`patchwork`](https://cran.r-project.org/web/packages/patchwork/index.html))
are available on CRAN. Others
([`TwoSampleMR`](https://mrcieu.github.io/TwoSampleMR/),
[`gsmr`](https://yanglab.westlake.edu.cn/software/gsmr/)) you will need
to install from source.

See the [Example](#example) below to get started.

## Usage

#### 1. Input files

To run the analysis you will need to include summary statistics for the
exposure and outcome of interest.

These should be in the form of a `data.frame`, and must contain the
following columns (case-sensitive): `rsid`, `effect_allele`,
`other_allele`, `beta`, `se`, `pval`.

`eaf` is not strictly required but this information is necessary to
perform MAF filtering and/or harmonisation of palindromic variants (see
below).

If your files contain additional columns, e.g. `chr`, `pos`, this is not
a problem.

#### 2. Preprocessing

##### Variant filtering

The function includes options for filtering variants in the exposure
data set.

1.  P-value thresholding `p_thresh =` : select variants that are
    associated with the exposure below a certain P-value threshold
    (default = 5e-8)

2.  MAF filtering `maf_thresh =`: select variants above a given minor
    allele frequency (MAF; default = 0.01)

3.  Selection of cis variants `cis =`: select cis variants within or
    around the gene region of interest (default = FALSE). You will need
    to specify the chromosome on which the gene is located, the starting
    and ending position of the gene, and define an optional flanking
    region using `gene_chr`, `gene_start`, `gene_end`, and `flank_kb`.

##### Harmonisation

The function will perform harmonisation of the exposure and outcome data
sets. You can read more about harmonisation in this
[`paper`](https://elifesciences.org/articles/34408).

Harmonisation is generally straightforward but requires an action for
dealing with palindromic SNPs – SNPs that have the same alleles on the
forward strand as the reverse strand (e.g. C/G on forward and G/C on
reverse). To read about the function options for dealing with
palindromic variants, please see the guidance provided
[here](PALINDROMES.md).

To harmonise the data, the function needs to match the variants from the
exposure and outcome data sets. By default, this is done based on the
rsid. However, if you wish to harmonise based on chromosome and
position, please see the guidance provided [here](CHRPOS.md).

##### Clumping

The function will optionally perform LD clumping of the data. You will
need to specify the distance for clumping (in kilobases; default=250),
the r-squared value for clumping (default = 0.01), and the reference
population (default = ‘EUR’) using the options `clump_kb`, `clump_r2`,
and `clump_pop`.

Note: By default, clumping is performed after harmonisation of the data,
which means variants will be excluded from clumping if they are not
found in the outcome data. If you wish to clump the data before
consideration of the outcome data, you can use `clump_only=T` to clump
the variants in the exposure data and output only a list of clumped
variants.

#### 3. Analysis

The pipeline will implement the following MR methods:

- [`Inverse variance weighted`](https://arxiv.org/abs/1512.04486)(fixed
  effects & random effects)
- [`Generalised inverse variance weighted`](https://onlinelibrary.wiley.com/doi/full/10.1002/gepi.21758)(fixed
  effects & random effects)
- [`Weighted median`](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4849733/)
- [`Weighted mode`](https://pubmed.ncbi.nlm.nih.gov/29040600/)
- [`MR-Egger`](https://academic.oup.com/ije/article/44/2/512/754653)
- [`MR-PRESSO`](https://www.nature.com/articles/s41588-018-0099-7)
- [`GSMR with HEIDI-outlier filtering`](https://www.nature.com/articles/s41467-017-02317-2)

You can read a high-level summary of these methods in this
[`paper`](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9123217/)

There are some specific considerations for some methods detailed below.

*Generalised IVW and GSMR*

> For Generalised IVW and GSMR, you will need to supply an LD matrix.
> For smaller numbers of SNPs (i.e. \<500), there is an option
> (`ld_mat = "try_gen"`) to request that the function create an LD
> matrix for you. However, there is no guarantee that all SNPs will be
> found in the reference panel (1000G V3).

*MR-PRESSO*

> MR-PRESSO recommends the number of bootstraps of the empirical
> distribution to be equal to nsnps/0.05 (e.g. 10,000 bootstraps for 500
> SNPs). The function will default to this recommended number, though a
> large number of bootstraps can cause significant increases in runtime.
> You can choose a smaller number of bootstraps, but the results may not
> be reliable.

## Example

##### Install packages

To get started, you will need to install the necessary packages:

``` r
# CRAN

install.packages("rio")
install.packages("tidyverse")
install.packages("data.table")
install.packages("MendelianRandomization")
install.packages("testthat")
install.packages("grid")
install.packages("gridExtra")
install.packages("ieugwasr")
install.packages("patchwork")
install.packages("remotes")
install.packages("survey")
install.packages("markdown")
install.packages("pacman")

# From source

remotes::install_github("MRCIEU/TwoSampleMR")
install.packages("https://yanglab.westlake.edu.cn/software/gsmr/static/gsmr_1.1.1.tar.gz",repos=NULL,type="source")
```

##### Load packages

After installing, you can load all packages using the
[`pacman`](https://cran.r-project.org/web/packages/pacman/index.html)
package:

``` r
if (!require(pacman)) install.packages("pacman")
p_load(rio, tidyverse,data.table, MendelianRandomization,
       testthat, grid,gridExtra,ieugwasr,patchwork,TwoSampleMR,gsmr )
```

##### Load function

Next, you will need to load the function and its dependencies:

``` r
# Main function
devtools::source_url("https://raw.githubusercontent.com/adamkvonende/MRpipeline/master/Scripts/perform_mr.R")
# Dependencies
devtools::source_url("https://raw.githubusercontent.com/adamkvonende/MRpipeline/master/Scripts/Functions/LiftSummaryStats.R")
devtools::source_url("https://raw.githubusercontent.com/adamkvonende/MRpipeline/master/Scripts/Functions/Map_rsid_from_chrpos.R")
```

##### Load exposure and outcome data

For this example, we will investigate the effects of Apolipoprotein C3
(ApoC3) on coronary heart disease (CHD) in Europeans. The first step is
to load the exposure and outcome data. For the purposes of illustration,
the exposure variants have been pre-filtered to P\<0.001.

``` r
# Load exposuredata: summary statistics for ApoC3

ss_apoc3 <- rio::import("https://zenodo.org/record/8090662/files/apoc3_summary_stats_pless001.csv?download=1")

# Load outcome data: summary statistics for CHD

ss_chd <- rio::import("https://zenodo.org/record/8090662/files/chd_summary_stats.csv?download=1") 
```

- **Example A**

First, we will an MR of ApoC3 on CHD using all variants significantly (P
\< 5e-8) associated with circulating ApoC3. We will perform clumping and
generate results for all methods. Since we need an LD matrix for some
methods, we will ask the function to generate one.

Note that we will run with fewer PRESSO bootstraps to save time.

``` r
results1 <- perform_mr(ss_exposure=ss_apoc3, # exposure data set
                       ss_outcome=ss_chd, # outcome data set
                       cis = F, # don't perform cis filtering (default)
                       maf_thresh=0.01, # filter on MAF <0.01 (default)
                       p_thresh=5*10^-8, # filter on P<5e-8 (default)
                       harmo_action=3, # drop palindromic variants
                       clump=T, # clump the variants (default)
                       clump_kb=250, # clump using 250k (default)
                       clump_r2 = 0.01, # clump r-squared = 0.01 (default)
                       clump_pop="EUR", # clump using European reference (default)
                       name_exp="Apoc3", # name the exposure variable
                       name_outcome="CHD", # name the outcome variable
                       which_methods="all", # perform all methods
                       ld_mat = "try_gen", # generate an LD matrix
                       ld_mat_pop="EUR", # create matrix using European reference
                       presso_nboot = 1000, # use 1000 bootstraps
                       which_plots="all", # plot all methods
                       outcome_binary=T) # specify that the outcome in binary (default)
```

Note that the same analysis can also be run with much less code as many
of the parameter values are defaults:

``` r
results1 <- perform_mr(ss_exposure=ss_apoc3, ss_outcome=ss_chd,name_exp="Apoc3", 
                       name_outcome="CHD", ld_mat = "try_gen", presso_nboot = 1000) 
```

<details>
<summary>
Show log for Example A
</summary>

    ## *** The function is now preprocessing the data for analysis ***
    ## 
    ##  5882 variants were retained in the exposure data after removing variants with p-value >= 5e-08
    ## 
    ##  5598 variants were retained in the exposure data after removing variants with maf < 0.01
    ## 
    ##  4684 variants were retained in the exposure data after removing variants with missing or identical alleles
    ## 
    ##  *** The function is now performing harmonisation ***
    ## 
    ##  65 variants in the exposure data set have missing rsids and will be dropped before harmonisation
    ## 
    ##  1486 SNPs were not found in the outcome data set according to rsid. Consider finding proxies. See 'harmo_exclude' for the list.
    ## 
    ##  There were 420 palindromic variants; these were dropped
    ## 
    ##  6 SNPs were excluded due to incompatible alleles. See 'harmo_exclude' for the list
    ## 
    ##  There are 2707 variants available for analysis after harmonisation
    ## 
    ##  *** The function is now performing clumping on the harmonised data***
    ## 
    ##  46 variants were retained after clumping at r2 = 0.01
    ## 
    ##  *** The function is now performing MR on the harmonised data ***
    ## 
    ##  The function is attempting to generate an LD matrix for the harmonised data...
    ## 
    ##  The LD matrix has been successfully created.
    ## 
    ##  Running MR-PRESSO using 1000 bootstraps specified by the user.
    ## 
    ##  *** The function has finished performing MR for Apoc3 and CHD ***

</details>

<br>

`perform_mr()` returns a named list containing the following objects:

- `results` – table of causal estimates for each method
- `het` – table of heterogeneity statistics
- `exclude` – list of outlier variants identified and method used to
  identify
- `harmo_exclude` – list of variants excluded in harmonisation and
  reason for exclusion
- `ld_mat` – LD matrix used in analysis
- `for_plot` – forest plot of results

Please note that the forest plot is provided as a rough visual aid and
is not designed to accommodate every use case. You can use the results
table provided to create customised plots.

You can access any of these objects using the ‘\$’ operator.

``` r
results1$results
results1$het
results1$exclude
results1$for_plot
```

<details>
<summary>
Show `results` table
</summary>
<p align="center">

<img src="https://i.ibb.co/6WCnPX6/results-table.png" alt="Description" style="width:70%;"/>

</p>
</details>
<details>
<summary>
Show `het` table
</summary>
<p align="center">

<img src="https://i.ibb.co/82zqFsg/het-table.png" alt="Description" style="width:40%;"/>

</p>
</details>
<details>
<summary>
Show `exclude` table
</summary>
<p align="center">

<img src="https://i.ibb.co/ckkYRNS/exclude-table.png" alt="Description" style="width:20%;"/>

</p>
</details>
<details>
<summary>
Show forest plot
</summary>
<p align="center">

<img src="https://i.ibb.co/C74ZcN9/forest-plot.png" alt="Description" style="width:70%;"/>

</p>
</details>

<br>

- **Example B**

Next we will run a cis-MR analysis of ApoC3 and CHD using all variants
significantly (p \< 5x10^-8) associated with ApoC3 concentration that
reside within 100kb of the APOC3 gene. We will again perform clumping,
but only generate results for ‘basic’ methods.

``` r
results2 <- perform_mr(ss_exposure=ss_apoc3, ss_outcome=ss_chd, 
                       cis = T, gene_chr=11, gene_start=116700422,gene_end=116703788,
                       flank_kb=100, clump=T, 
                       name_exp="Apoc3", name_outcome="CHD",
                       which_methods="basic")
```

<details>
<summary>
Show log for Example B
</summary>

    ## *** The function is now preprocessing the data for analysis ***
    ## 
    ##  5882 variants were retained in the exposure data after removing variants with p-value >= 5e-08
    ## 
    ##  5598 variants were retained in the exposure data after removing variants with maf < 0.01
    ## 
    ##  4684 variants were retained in the exposure data after removing variants with missing or identical alleles
    ## 
    ##  174 cis variants were found in the exposure data and carried forward for analysis
    ## 
    ##  *** The function is now performing harmonisation ***
    ## 
    ##  16 SNPs were not found in the outcome data set according to rsid. Consider finding proxies. See 'harmo_exclude' for the list.
    ## 
    ##  There were 21 palindromic variants; these were dropped
    ## 
    ##  
    ## There are 137 variants available for analysis after harmonisation
    ## 
    ##  *** The function is now performing MR on the harmonised data ***
    ## 
    ##  Apoc3 - CHD
    ##  *** The function has finished performing MR for Apoc3 and CHD ***

</details>
<details>
<summary>
Show `results` table
</summary>
<p align="center">

<img src="https://i.ibb.co/9nJWTJN/results-table2.png" alt="Description" style="width:70%;"/>

</p>
</details>
<details>
<summary>
Show forest plot
</summary>
<p align="center">

<img src="https://i.ibb.co/VS6km5Z/forest-plot2.png" alt="Description" style="width:70%;"/>

</p>

## Contact

If you have any questions or suggestions please don’t hesitate to
contact [Adam Von Ende](mailto:adam.vonende@ndph.ox.ac.uk) or [Elsa
Valdes-Marquez](mailto:elsa.valdes.marquez.vonende@ndph.ox.ac.uk).
