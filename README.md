# spruceGBS

Please note, this is a work in progress.

---

# Background

This repository contains workflows and analyses for a large-scale genotyping-by-sequencing (GBS) study in _Picea_ (spruce). The dataset includes ~1,500 individuals sampled from ~150 locations spanning the Atlantic to Pacific coasts and from the Arctic Ocean to the Adriatic Sea, covering _Picea abies_ (Norway spruce), _P. obovata_ (Siberian spruce), and possibly _P. koraiensis_. These are long-lived, outcrossing forest trees with large, repetitive genomes (~20 Gb, > 70% TEs) and an incomplete reference assembly, factors that complicate short-read variant discovery.

---

# Motivation 

We initially applied a standard, literature-based variant calling pipeline (`fastQC` → `BWA-MEM` → indel realignment with legacy `GATK` → `samtools mpileup` → basic SNP filters). While this approach recovered the expected broad-scale population structure, it also suggested a surprising spatial distribution of genetic diversity (e.g.,  extremely low diversity in some regions with implausibly high (F<sub>is</sub> > 0.40) inbreeding estimates), suggesting the generic pipeline was inadequate for this dataset.

This repository develops and evaluates a pipeline tailored to the specific challenges of the spruce GBS data, including:

- multiple sequencing libraries (some HiSeq, some NovaSeq) with varying DNA quality and input quantity;
- partial confounding between geography and sequencing library;
- a highly fragmented (10 million scaffolds), incomplete (~ 60%) reference assembly (circa 2013);
- potential reference bias, paralog collapse, null alleles, and allelic dropout across divergent lineages and hybrid zones.

---
# Objectives

The goal is to produce more reliable estimates of genetic diversity and structure by systematically testing parameter choices and filters, with an emphasis on transparency and reproducibility. I have tried to write didactic READMEs, but this, in particular, is a work in progress. 

Additionally, I am using this project to explore ways to make the implementation of customized pipelines more accessible. This includes things like the statistical summaries used in the parameter sweep to make interpretation more scalable and less qualitative, the interpretable machine learning approach using gradient boosted decision trees for variant filtration, and scripts to make parsing the results from `ANGSD` in particular less angst-inducing. The most developed of these tools are:
- [`mqforensics`](https://github.com/lxsllvn/spruceGBS/tree/main/mqforensics)  - a C program for characterizing mapping quality distributions, including capped (`-C`) qualities, INDELs, clipped bases, mismatches and more from BAM alignments.
- ['beagle-utils'](https://github.com/lxsllvn/spruceGBS/tree/main/beagle-utils) - a Python package for manipulating and subsetting read count matrices and beagle-formatted genotype likelihoods produced by `ANGSD`.

---

# A note on the reference assembly

Most of the pipeline should be applicable to other datasets, but handling the size of the *P. abies* *v*. 1.0 assembly requires some case-specific steps that would be unnecessary with, for example, a mammal. 

Reads were mapped to the full assembly to maximize the chance of detecting blah blah blah. Then, we progressively reduced the number of analyzed scaffolds, first during the [reduced reference preparation](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref) to a size manageable by [GATK's IndelRealigner](https://github.com/lxsllvn/spruceGBS/tree/main/04_realignment) (unfortunately deprecated in v. 3.6) and more aggressively during the [ANGSD site discovery](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery). The [ANGSD parameter sweep](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep) was conducted over 100 Mbp from the longest scaffolds. 

Creating references for the ANGSD analyses was a three-step process (which we repeat for each of the three domains, or major *P. abies*-*obovata* lineages, covered by the dataset). For a smaller and/or more contiguous assembly, creating a special [experimental reference](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#scaffold-selection) for the parameter sweep may be unnecessary and [site discovery](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery) could be carried out over the entire reference at once, rather than on 23 reference subsets. 

---

# Pipeline contents 

# [Step 1: initial pre-processing and alignment](https://github.com/lxsllvn/spruceGBS/tree/main/01_read_alignment)
  * [Read pre-processing and quality control](https://github.com/lxsllvn/spruceGBS/tree/main/01_read_alignment#read-pre-processing-and-quality-control)
  * [Alignment](https://github.com/lxsllvn/spruceGBS/tree/main/01_read_alignment#pre-processing-and-quality-control)

# [Step 2: reduced reference preparation](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref)
  * [Scaffold search](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref#scaffold-search)
  * [Reduced reference preparation](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref#reduced-reference-preparation)
  * [Identify target regions for analysis](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref#identify-target-regions-for-analysis)
  * [BAM intersections](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref#bam-intersections)
  
# [Step 3: initial sample quality control](https://github.com/lxsllvn/spruceGBS/tree/main/03_initial_qc)
  * [Sequencing depth and breadth per sample](https://github.com/lxsllvn/spruceGBS/tree/main/03_initial_qc#sequencing-depth-and-breadth-per-sample)
  * [Visualization and removal](https://github.com/lxsllvn/spruceGBS/tree/main/03_initial_qc#visualization-and-removal)
  
# [Step 4: indel realignment](https://github.com/lxsllvn/spruceGBS/tree/main/04_realignment)
  * [Important notes](https://github.com/lxsllvn/spruceGBS/edit/main/04_realignment/README.md#important-notes)
  * [Find targets for indel realignment](https://github.com/lxsllvn/spruceGBS/edit/main/04_realignment/README.md#find-targets-for-indel-realignment)
  * [Indel realignment](https://github.com/lxsllvn/spruceGBS/edit/main/04_realignment/README.md#indel-realignment)
  
# [Step 5: ANGSD parameter sweep](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep)
  * [Objectives](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#objectives)
  * [Filter parameter discussion](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#filter-parameters)
  * [Experimental design](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#experimental-design)
    * [Sample selection](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#sample-selection)
    * [Scaffold selection](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#scaffold-selection)
  * [Site and population-level statistics](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#site-and-population-level-statistics)
  * [PCA and dd-RDA](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#rda-on-principal-components)
  * [Individual heterozygosity](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#individual-heterozygosity)
    
# [Step 6: analyze parameter sweep](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results) 
  * [Summary](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#summary)
  * [Statistical methods](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#statistical-methods)
  * [Core parameter combination results](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#core-parameter-combinations)
    * [Individual heterozygosity](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#individual-heterozygosity)
      * [Mixed-effect models](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#mixed-effect-models)
      * [Within-population variation](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#within-population-variation)
    * [RDA and PCA](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#rda-and-pca)
  * [Extended parameter results](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#extended-parameter-combinations)
    * [RDA and PCA](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#rda-and-pca-1)
    * [Individual heterozygosity](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#individual-heterozygosity-1)
      * [Mixed-effect models](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#mixed-effect-models-1)
      * [Within-population variation](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#within-population-variation-1)
    * [Genetic diversity results](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#genetic-diversity)
   * [Scripts](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#scripts)
   * [Function index](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#analysis_functionsr-index)

# [Step 7: find and filter sites](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery)
* [Create ANGSD reference assemblies](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#create-angsd-reference-assemblies)
* [Calculate genotype likelihoods](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#create-angsd-reference-assemblies)
* [Site filtering](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#site-filtering)
  * [Interpretable machine learning](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#interpretable-machine-learning)
  * [Feature engineering](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#feature-engineering)
