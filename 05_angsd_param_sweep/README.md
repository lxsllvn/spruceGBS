# Overview

Implements **Step 5: ANGSD parameter sweep** of the spruceGBS pipeline.

---

# Contents

* [Objectives](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#objectives)
* [Filter parameter discussion](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#filter-parameters)
* [Experimental design](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#experimental-design)
* * [Sample selection](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#sample-selection)
  * [Scaffold selection](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#scaffold-selection)
* [Parameter sweep: site and population-level statistics](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#site-and-population-level-statistics)
* [Parameter sweep: PCA and MANOVA](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#manova-on-principal-coordinates)
* [Parameter sweep: individual heterozygosity](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#individual-heterozygosity)

---

# Scripts

* **`param_exp_cap_mapq.py`**:
* **`param_exp_sample_selection.R`**: 
* **`param_exp_ref_step1.sh`**: 
* **`param_exp_ref_step2.sh`**:
* **`param_exp_ref_step3.sh`**:
* **`param_exp_popstats.sh`**:
* **`param_exp_summarize_popstats.R`**:
* **`param_exp_PCA.sh`**:
* **`param_exp_manova.R`**:
* **`param_exp_indvhet.sh`**: 
---

# Objectives

Variant filters have a profound impact on downstream estimates of population structure and diversity, with the potential to either mitigate or exacerbate biases introduced during wet-lab protocols. Choosing optimal filter parameters is particularly critical for this dataset, which presents several challenges: a non-model organism with a large and repetitive genome, multiple diverged lineages, sequencing across different libraries and platforms, and alignment to a highly fragmented reference genome.

Our objective was to assess the influence of filter parameter settings across three key analyses, aiming to minimize technical artifacts (e.g., batch effects, allelic dropout) while preserving known biological signal:

1. **Population-level site statistics:**
    - Evaluated the impact of parameter choices on H<sub>e</sub>, H<sub>o</sub>, F, π, Θ<sub>w</sub>, and MAF 
    - Identified settings that led to biologically implausible estimates, if any.
      
2. **PCA and MANOVA on principal coordinates:**
    - Tested the variance explained by sequencing library and geographic region.
      
3. **Individual heterozygosity estimates:**
    - Assessed parameter effects on individual variation in mixed-library populations.

We conducted the analyses over a grid of parameter settings:

- **Base Quality (`-minQ`)**: 20, 30
- **Mapping Quality (`-minMapQ`)**: 20, 30, 40
- **BAQ Model (`-baq`)**: 0 (off), 1 (simple), 2 (extended)
- **Mapping Quality Capping Coefficient (`-C`)**: 0, 50
- **Per-Library Call Rate Filters**: 40%, 50%, and 60%

To further refine parameter choices, we extended analyses to a wider range of `-C` values after identifying promising combinations.

Below, we discuss the potential effects of each filter parameter. For implementation details, skip to [Experimental-design](https://github.com/lxsllvn/spruceGBS/edit/main/05_angsd_param_sweep/README.md#experimental-design).

---

# Filter parameters 

## Minimum base quality (`-minQ`)

Our dataset includes libraries sequenced on the HiSeq and NovaSeq. We attempted to mitigate platform differences during read quality control, but identifying an appropriate minimum base quality must account for three important differences: 
1. **2- vs. 4-channel chemistry**. Unlike the Hiseq, where each base has its own dye, NovaSeq bases are encoded by combinations of two dyes, and low-quality or missing bases are often called as a G.
2. **Binned vs. continuous Q‐scores**. The HiSeq reports continuous Phred scores, whereas the Novoseq collapses Q scores in two eight discrete bins. As a result, the same reported Q score does not necessarily indicate the same error probability. 
3. **Underlying error spectra**. Platforms differ in both the type of sequencing errors they produce and the context of these errors. Hiseq reads show pronounced cycle-dependent quality drop-offs and errors tend to cluster in just a few sequence contexts (e.g. [Minoche et al. 2011](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-11-r112)). NovaSeq reads are enriched for poly-g artifacts as a result of their 2-channel chemistry, but other error rates tend to be more uniform across sequence contexts (e.g. [Ma et al. 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1659-6)).

In addition, samples differ in the quality and amount of their starting DNA. Degraded DNA may show more deamination artefacts and lower per‐base quality, and the minimum Q must balance false positives against the risk of exacerbating library-specific missingness. 

## Minimum mapping quality (`-minMapQ`)

The *Picea* nuclear genome is > 70% TEs, and while the GBS protocol enriches the non-repetitive fraction of the genome, this presents a problem for ensuring that reads mapping to the same coordinates actually originated from the same place in the genome. Our dataset spans at least diverged populations, whereas the reference assembly originates from a single tree in central Sweden. Reads from diverged populations can have systematically lower MQs, even if they are correctly mapped to their homologs in the reference. The choice of minimum mapping qualities must therefore balance against spurious alleles without systematically reducing coverage and biasing allele frequencies in non-reference populations. 

## Base alignment quality (`-baq`) model

Base Alignment Quality (BAQ) uses a small, banded dynamic-programming realignment around each read’s original CIGAR string to compute per-base posterior probabilities of misalignment (Phred-scaled). BAQ annotates existing alignments by adjusting the QUAL field or writing a BQ/ZQ tag, but it does not rewrite the CIGAR or move bases.

ANGSD has three BAQ options. BAQ = 0 turns it off, but the [documentation](https://www.popgen.dk/angsd/index.php/Input#Arguments) is [unclear](
https://github.com/ANGSD/angsd/issues/97#issuecomment-327461936) on what BAQ = 1 and BAQ = 2 actually do. The confusion probably comes from `samtools` using 'partial' to refer to heuristics for selecting *which* reads to realign and using 'extended' in older versions to refer to *how* the qualities were adjusted. 

In any event, the source code (annotations and any errors below are mine) shows that:
* `-baq 1` ('partial' according to ANGSD but maybe 'simple' or 'local' are less misleading?): qualities for bases within each contiguous match/mismatch block (CIGAR M/=/X) that are flagged by the DP as an insertion, deletion, or shifted alignment are zeroed; the others are capped by the DP’s posterior Phred score.
* `-baq 2` (extended): does the same zero/cap pass, then for each match block computes left-and-right running maxima over those provisional scores and takes the minimum at each position. This smooths low-confidence spikes into their neighbors, boosting sensitivity to SNPs near indels.

BAQ is expected to reduce false positive SNPs around indels ([Li 2011](https://doi.org/10.1093/bioinformatics/btr076](https://doi.org/10.1093/bioinformatics/btr076))). In practice, different BAQ models may lead to markedly different site frequency spectrum, but whether this is also a more accurate site frequency spectrum may be unclear and potentially dependent on factors such as sequencing depth (see discussions [here](https://www.biostars.org/p/9466154/), [here](https://github.com/ANGSD/angsd/issues/97) and [here](https://github.com/samtools/bcftools/pull/1363#issuecomment-80204203).) 

```C
int sam_prob_realn_core(bam1_t *b, const char *ref, int flag) {
    // Extract core alignment and quality information
    bam1_core_t *c = &b->core;
    uint8_t *qual = bam1_qual(b);  // Original base QUAL array

    // Decode BAQ flags from 'flag'
    int apply_baq  = flag & 1;        // bit 0: apply BAQ (write ZQ or BQ)
    int extend_baq = (flag >> 1) & 1; // bit 1: use extended (smoothed) BAQ

    // 1) Skip unmapped / zero-length reads
    if ((c->flag & BAM_FUNMAP) || c->l_qseq == 0)
        return -1;

    // 2) Handle existing BQ/ZQ tags (convert, skip, or exit early)
    //    (omitted here for brevity)

    // 3) Compute reference window [xb, xe) covering the read
    //    by walking the CIGAR, also compute DP bandwidth 'bw'
    //    to limit the DP matrix size. (code omitted)

    // 4) Allocate and fill arrays for DP:
    //    s[]    = read bases (encoded 0-3),
    //    r[]    = reference bases over [xb, xe),
    //    state[]= DP-chosen state at each read base,
    //    q[]    = DP posterior Phred quality per base,
    //    bq[]   = final BAQ adjustments (to write tag)
    uint8_t *s     = calloc(c->l_qseq, 1);
    uint8_t *r     = calloc(xe - xb, 1);
    uint8_t *bq    = calloc(c->l_qseq, 1);
    int     *state = calloc(c->l_qseq, sizeof(int));
    uint8_t *q     = calloc(c->l_qseq, 1);

    // 5) Run banded "glocal" DP (posterior realignment):
    //    - Smith–Waterman style DP under an HMM modeling matches,
    //      mismatches, insertions, and deletions.
    //    - Fills 'state[i]' with low 2 bits = (0:match,1:ins,2:del,3:mm),
    //      high bits = reference offset along the DP diagonal.
    //    - Computes 'q[i]' as Phred-scaled posterior confidence that
    //      read-base i is correctly aligned to the ref.
    kpa_par_t conf = kpa_par_def;
    kpa_glocal(r, xe - xb, s, c->l_qseq, qual, &conf, state, q);

    // 6) Partial vs. Extended BAQ based on 'extend_baq' flag
    if (!extend_baq) {
        //    Partial BAQ:
        //    For each match block in the CIGAR (M,=,X ops):
        //      - If DP state[i] != match, zero bq[i] (ignore base)
        //      - Else bq[i] = min(original cap, q[i])
        for (int k = 0, x = c->pos, y = 0; k < c->n_cigar; ++k) {
            int op = cigar[k] & 0xf, l = cigar[k] >> 4;
            if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                for (int i = y; i < y + l; ++i) {
                    // Check if DP aligns base i to the expected ref pos
                    if ((state[i] & 3) != 0 ||
                        state[i] >> 2 != (x - xb) + (i - y)) {
                        bq[i] = 0;
                    } else {
                        // Cap by DP posterior
                        bq[i] = q[i] < bq[i] ? q[i] : bq[i];
                    }
                }
                x += l;
                y += l;
            }
            // handle I/S/D ops advancing x or y (omitted)
        }
        // Convert to ASCII offset for BQ/ZQ tag
        for (int i = 0; i < c->l_qseq; ++i)
            bq[i] = qual[i] - bq[i] + 64;

    } else {
        //    Extended BAQ:
        //    1) Perform same zero/cap logic into bq[] as partial BAQ.
        //    2) Compute left-to-right and right-to-left running maxima
        //       over provisional bq[], then bq[i] = min(left[i], right[i]).
        uint8_t *left = calloc(c->l_qseq, 1);
        uint8_t *rght = calloc(c->l_qseq, 1);
        // ...compute left[] and rght[] sweeps based on provisional bq[]...
        for (int i = 0; i < c->l_qseq; ++i)
            bq[i] = left[i] < rght[i] ? left[i] : rght[i];
        // Apply ASCII offset
        for (int i = 0; i < c->l_qseq; ++i)
            bq[i] += 64;
        free(left);
        free(rght);
    }

    // 7) Write out BQ or ZQ tag:
    //    - If apply_baq: subtract (bq[i]-64) from QUAL[i] and append ZQ tag.
    //    - Else: leave QUAL intact and append BQ tag.
    //    (Tag-writing code omitted)
}
```

## Mapping quality capping coefficient (`-C`)

High mapping qualities do not necessarily indicate good alignments.  A read may align uniquely and earn MQ = 60 even if it carries dozens of high-quality mismatches (e.g. 20 × Q ≥ 20 → ΣQ ≥ 400) or large soft-clips (e.g. 50 × Q≈20 → ΣQ≈1,000).

The `-C` parameter is used to downgrade mapping qualities in cases like these. When invoked, the capped MQ function (`sam_cap_mapq` in the source code) first computes an error penalty *T* from the sum of high-quality mismatches (SubQ) and clipped-base penalties (ClipQ).  

Then, if *T* ≤ *C*, the new MQ is $= C×sqrt((C−T)/C)$ up to a maximum of *C*  and MQ = 0 if *T* > *C*. 

This down-weights moderately noisy alignments and eliminate the worst ones. Choosing *C* below the read aligner’s maximum MQ (e.g. < 60) will increasingly penalize true heterozygote reads (which carry one genuine mismatch), whereas higher values of  *C*  filter only the most egregious misalignments. 

`param_exp_cap_mq.py` illustrates the effects of adjusting MQs. Given a bam and a `-C` value, it calculates the new MQ and reports the original MQ, the read name and length, SubQ, ClipQ, *T*, and AS1, the primary alignment score. 

## Library call rates

Filtering on overall site call-rate can hide, or even worsen, systematic missingness that tracks sequencing library rather than biology. When one library routinely fails to call certain loci (due to GBS fragment size-selection, DNA quality, GC/AT bias, or simply lower depth), those sites carry a “library fingerprint” of missing data. Downstream tools (PCA, ADMIXTURE) that impute or mean-replace missing genotypes will then recover batch effects instead of true population structure.

---

# Experimental design

## Sample selection

We identified samples from each domain where library effects should be most distinguishable from real genetic structure. Specifically, these samples were:

1.  ≥ 25% percentile in mapped reads and scaffold coverage to exclude potentially unrepresentative marginal samples
2. from populations (i.e., a single collection site) with ≥ 5 sampled trees, 

and were:

3. from stands where ≤75% of trees were sequenced in a single library (56% on average), or
4. within 100 km (≤ 250 km for Siberia) of a population with a different majority library.

Given the outcrossing mating system of *Picea*, trees from the same collection site are unlikely to deviate from Hardy-Weinberg equilibrium or show strong variation in individual heterozygosity. Genetic structure develops slowly in *Picea*, making batch effects a likely source of any apparent differentiation between stands in relatively close proximity. 

### `param_exp_sample_selection.R` usage
```R

# Step 1: Filter and annotate samples
samples <- load_and_filter_samples(path = "intersected_bam_mapping_summary.tsv", 
                                   meta_path = "sequenced_samples_metadata.csv", 
                                   read_quantile = 0.25, scaff_quantile = 0.25)

# Step 2: Identify mixed-library pops 
mixed_pops <- summarize_populations(samples = samples, 
                                    min_samples = 5, 
                                    max_prop = 0.8,
                                    mixed_only=TRUE)

# Step 3: Build population lookup
lib_lookup <- summarize_populations(samples = samples, 
                                    min_samples = 5, 
                                    mixed_only = FALSE, 
                                    keep_cols = "domain")

# Step 4: Find population pairs using geodesic distance matrix
pairs <- find_population_pairs(geo_matrix_path = "pops_geodesic_dist.txt",
                               lib_lookup = lib_lookup)

# Step 5: Build final status table
status_table <- build_population_status(mixed_tbl = mixed_pops, 
                                        paired_tbl = pairs,
                                        lib_lookup = lib_lookup)
```
**Inputs**
* `intersected_bam_mapping_summary.tsv`: tab-delimited summary of mapped reads and scaffold coverage from [initial QC](https://github.com/lxsllvn/spruceGBS/tree/main/03_initial_qc)
* `sequenced_samples_metadata.csv`: sample metadata table; must include: `bam_code`, `pop_code`, `domain`, `region`, `library`, `latitude`, `longitude` 
* `pops_geodesic_dist.txt`: tab-delimited distance matrix of pairwise geodesic (km) distances between populations

**Outputs**
* `angsd_parameter_exp_pops.csv`: summary table of populations used in the parameter sweep experiment 

## Scaffold selection

Finding sites that yield enough SNPs for analysis while minimizing resources during the parameter sweep requires some setup. First, I took 100 Mbp from the longest scaffolds, divided this into 10 reference subsets (`param_exp_ref_step1.sh`). Then, I ran ANGSD on each reference subset to find sites with < 50% missing data, using samples from each domain separately (`param_exp_ref_step2.sh`). Lastly, I created the experimental references from these sites, resulting in an indexed fasta and the ANGSD site and region files for each domain, e.g. siberia_experiment_ref.fa, siberia_experiment_ref_regions, and siberia_experiment_ref_sites (`param_exp_ref_step3.sh`). 

### Create reference subsets

select 100 Mbp from the longest putatively non-contaminant scaffolds (i.e., those without "putative contaminant" in their fasta headers), divided this into 10 subsets, and created a indexed reference genome and ANGSD region and site files

#### `param_exp_ref_step1.sh` usage

```
#!/bin/bash
sbatch param_exp_ref_step1.sh \
    "$SPRUCE_PROJECT/parameter_testing/exp_ref" \
    "$SPRUCE_PROJECT/ref/picea_newref_target_regions.bed" \
    "$SPRUCE_PROJECT/ref/ref_putative_contaminants.bed" \
    "$SPRUCE_PROJECT/ref/picea_newref.fa"
```

**Inputs**
* `$1` – output directory for the reference subsets
* `$2` – bed file of [target regions](https://github.com/lxsllvn/spruceGBS/blob/main/02_reduced_ref/)
* `$3` – bed file of potential contaminant scaffolds; parsed from `Pabies1.0-genome.fa` fasta headers
* `$4` – path to [reduced reference genome](https://github.com/lxsllvn/spruceGBS/blob/main/02_reduced_ref/)

**Outputs**
* `experiment_ref_pt_a{a..j}.fa` - reference subset fasta
* `experiment_ref_pt_a{a..j}.fai` - fasta index
* `target_scaff_pt_a{a..j}_regions` - ANGSD region file
* `target_scaff_pt_a{a..j}_sites` - ANGSD sites file

### Domain-level site discovery

For each domain, I then calculated the read count matrix of each reference subset using the most lenient set of filter combinations (`-baq 0 -C 0 -minQ 20 -minMapQ 20`) and required a site to be present in at least half the samples.

#### `param_exp_ref_step2.sh` usage

 Allocate ca. 12 gb memory for southern and siberia; ca. 24 for northern 

```
#!/bin/bash
# Loop through parameter sets a to j and submit each job to SLURM
for i in {a..j}; do
    sbatch "$SCRIPTS/angsd_param_exp_discovery.sh" \
        "$i" \
        "southern" \
        "$SPRUCE_PROJECT/parameter_testing/parameter_exp_southern_bamlist" \
        "$SPRUCE_PROJECT/parameter_testing/ex_ref/southern_site_discovery"
done
```

**Inputs**
* `$1` – the chunk index (e.g., `a`, `b`, ...)
* `$2` – the domain or region (e.g., `southern`)
* `$3` – path to the BAM list file (list of input BAMs)
* `$4` – output directory for results

**Outputs**
* `${DOMAIN}_target_scaff_pt_a{a..j}.counts.gz` - read count matrices for each chunk 
* `${DOMAIN}_target_scaff_pt_a{a..j}.pos.gz` - scaffold and position names for each chunk

**NB!** The `*.counts.gz` matrix does not have informative row or column names. Most, but not all, ANGSD output follows this convention. Some files with positional data from ANGSD use `chr` and `pos`, like `.pos.gz`, but others use `Chromo` and `position` and there's probably other versions as well. 

### Experimental reference preparation

Then, concatenate the pos.gzs, parse into a bed file and ANGSD region and site files, and build and index experimental reference  

#### `param_exp_ref_step3.sh` usage

```
#!/bin/bash
sbatch ${SCRIPTS}/angsd_param_exp_domain_refprep.sh \
	"southern" \
	"${SPRUCE_PROJECT}/parameter_testing/exp_ref/southern_site_discovery" \
	"${SPRUCE_PROJECT}/parameter_testing/exp_ref"
```

**Inputs**
- `$1` – the domain or region (e.g., `southern`)
- `$2` – directory containing `${DOMAIN}_target_scaff_pt_a{a..j}.counts.gz` and `${DOMAIN}_target_scaff_pt_a{a..j}.pos.gz` for each chunk; created by `param_exp_ref_step2.sh`
- `$3` – output directory for the reference files

**Output**
* An experimental reference comprising sites with < 50% missing data and associated files per domain:
  * `${DOMAIN}_experiment_ref.fa`
  * `${DOMAIN}_experiment_ref.fa.fai`
  * `${DOMAIN}_experiment_ref.bed`
  * `${DOMAIN}_experiment_ref_sites`
  * `${DOMAIN}_experiment_ref_regions`

---
  
# Site and population-level statistics

Uses the [experimental samples](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#sample-selection) and [experimental reference genome](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#scaffold-selection)

Calculates domain-level DP matrix and MAFs
Calculates H<sub>e</sub>, H<sub>o</sub>, F, π, Θ<sub>w</sub>, and MAF by locus per population. 

Call threshold and MAF filters are implemented in the R script, which uses the domain-level read depth matrix and MAFs. R script also prepares the data for the PCA and MANOVA.

## **`param_exp_popstats.sh`** usage

```bash
#!/bin/bash
sbatch $SCRIPTS/angsd_param_exp_sweep.sh \
"$SPRUCE_PROJECT/parameter_testing/exp_ref/southern_experiment_ref.fa" \
southern \
parameter_exp_southern_bamlist \
"${SPRUCE_PROJECT}/parameter_testing/exp_ref/southern_experiment_ref_regions" \
"${SPRUCE_PROJECT}/parameter_testing/exp_ref/southern_experiment_ref_sites" \
southern_results
```

**Inputs**
* `$1` – path to the [experimental reference genome](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#scaffold-selection) for a domain
* `$2` – domain name (e.g., southern)
* `$3` – list of BAM files in the domain (e.g., parameter_exp_southern_bamlist)
* `$4` – region file for the [experimental reference genome](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#scaffold-selection)
* `$5` – sites file for the [experimental reference genome](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#scaffold-selection) 
* `$6` – base output directory name 

Per-population BAM lists should be stored in "${DOMAIN}_populations/*.txt"

**Outputs** 
```bash
├── ${OUTDIR}/${DOMAIN}_${PARAM_ID}/  
├── ${POP}/${POP}/            # one folder per population 
│   ├── gl.arg                # run arguements
│   ├── gl.beagle.gz          # genotype likelihoods, beagle format
│   ├── gl.hwe.gz             # per-site inbreeding coefficients 
│   ├── gl.mafs.gz            # minor allele frequencies 
│   ├── gl.safs.gz            # binary site allele frequency likelihood
│   ├── gl.safs.idx           # binary postion index
│   ├── gl.safs.pos.gz        # binary scaffold names and positions
│   ├── gl.sfs                # site frequency spectrum
│   ├── saf2theta.thetas.gz   # binary per-site Θ estimates
│   ├── saf2theta.thetas.idx  # binary index file
│   ├── thetas.persite.txt    # neutrality test statistics π and Θ estimators
│   ├── thetas.summary.pestPG # per-scaffold summary π and Θ estimators
└── domain_sfs.counts.gz      # read count matrix, nrow sites X ncol samples
└── domain_sfs.pos.gz         # scaffold names and postions
└── domain_sfs.mafs.gz        # domain-level minor allele frequences
└── domain_sfs.beagle.gz      # domain-level genotype likelihoods in 
                              # beagle format, used for PCAngsd later
```

Results from each parameter combination are written to `${DOMAIN}_${PARAM_ID}`. The `${PARAM_ID}` name specifies the `-baq` model, `-C` coefficent, `-minQ` and `-minMapQ` used in the run; for example `baq0_C0_q20_mq20`, `baq0_C0_q20_mq30`, `baq0_C0_q20_mq40`, ... . Domain-level results live in this folder and population-level results each live in their own subdirectory. 

Once `param_exp_popstats.sh` is completed, `param_exp_popstats.R` implements the library call-rate filter, creates genotype likelihood files for each call rate for later use, and summarizes the results in two files: `$DOMAIN_angsd_param_summaries.csv` for all loci and `$DOMAIN_maf05_angsd_param_summaries.csv` for loci with a domain-level MAF > 0.5. 

## **`param_exp_popstats.R`** usage

```bash
Rscript "${SCRIPTS}/angsd_param_exp_summary.R" \
    "${SPRUCE_PROJECT}/parameter_testing/southern_results" \
    "${SPRUCE_PROJECT}/parameter_testing/sequenced_samples_metadata.csv" \
    "${SPRUCE_PROJECT}/parameter_testing/parameter_exp_southern_bamlist" \
    "southern"
```

**Inputs**
* `$1` – path to top-level results folder 
* `$2` – sample metadata table; must include: `bam_code`, `pop_code`, `domain`, `region`, `library`, `latitude`, `longitude`
* `$3` – list of BAM files in the domain (e.g., parameter_exp_southern_bamlist)
* `$4` – domain name (e.g., southern)


**Outputs**
```bash
├── ${OUTDIR}/${DOMAIN}_${PARAM_ID}/
│   ├── ${DOMAIN}_${PARAM_ID}_ct4.beagle.gz  # beagle genotype likelihoods for sites with library call rate > 40%
│   ├── ${DOMAIN}_${PARAM_ID}_ct5.beagle.gz  # library call rate > 50%
│   ├── ${DOMAIN}_${PARAM_ID}_ct6.beagle.gz  # library call rate > 60%
└── $DOMAIN_angsd_param_summaries.csv        # collected results for all loci
└── $DOMAIN_maf05_angsd_param_summaries.csv  # collected results for domain-level MAF > 0.05 loci
```

`${DOMAIN}_${PARAM_ID}_ct4.beagle.gz` are the input files for PCAngsd in the next step. 

---

# MANOVA on principal coordinates 

---

# Individual heterozygosity

---

# Dependencies

List required modules, software, or packages:

* 

---

