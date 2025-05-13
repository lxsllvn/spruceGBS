## Overview

Briefly describe the purpose of this directory (one or two sentences).

---

## Contents

* **`<script1>.sh`**: Short description of what this script does.
* **`<script2>.R`**: Short description of what this analysis or step does.
* **`...`**

---
# Objectives

Identify settings that minimize technical artifacts (e.g., batch effects, allelic dropout) while retaining true biological signal

- **Base quality**: `-minQ 20`, `30` 
- **Mapping quality**: `-minMapQ 20`, `30`, `40`
- **BAQ model**: `-baq 0`, `1`, `2`
- **Coefficient for capping mapping qualities**: `-C 0`, `50`
- **Library call rates** Per-library missing data filters (\< 40, 50, and 60%),

## Minimum base quality (`-minQ`)

Our dataset includes libraries sequenced on the HiSeq and NovaSeq. We attempted to mitigate platform differences during read quality control, but identifying an appropriate minimum base quality must account for three important differences: 
1.   **2- vs. 4-channel chemistry**. Unlike the Hiseq, where each base has its own dye, NovaSeq bases are encoded by combinations of two dyes, and low-quality or missing bases are often called as a G.
2.  **Binned vs. continuous Q‐scores**. The HiSeq reports continuous Phred scores, whereas the Novoseq collapses Q scores in two eight discrete bins. As a result, the same reported Q score does not necessarily indicate the same error probability. 
3. **Underlying error spectra**. Platforms can differ in the type of sequencing errors they produce (e.g.) and the context of these errors. Hiseq reads show pronounced cycle-dependent quality drop-offs and errors tend to cluster in just a few sequence contexts (e.g. [Minoche et al. 2011](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-11-r112)). NovaSeq reads are enriched for poly-g artifacts as a result of their 2-channel chemistry, but other error rates tend to be more uniform across sequence contexts (e.g. [Ma et al. 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1659-6)).

In addition, samples differ in the quality and amount of their starting DNA. Degraded DNA may show more deamination artefacts and lower per‐base quality, and the minimum Q must balance  false positives against the risk of exacerbating library-specific missingness. 

## Minimum mapping quality (`-minMapQ`)

The *Picea* nuclear genome is > 70% TEs, and while the GBS protocol enriches the non-repetitive fraction of the genome, this presents a problem for ensuring that reads mapping to the same coordinates actually originated from the same place in the genome. Our dataset spans at least diverged populations, whereas the reference assembly originates from a single tree in central Sweden. Reads from diverged populations can have systematically lower MQs, even if they are correctly mapped to their homologs in the reference. The choice of minimum mapping qualities must therefore balance against spurious alleles without systematically reducing coverage and biasing allele frequencies in non-reference populations. 

## Base alignment quality (`- baq`) model

Base Alignment Quality (BAQ) uses a small, banded dynamic-programming realignment around each read’s original CIGAR string to compute per-base posterior probabilities of misalignment (Phred-scaled). BAQ annotates existing alignments by adjusting the QUAL field or writing a BQ/ZQ tag, but it does not rewrite the CIGAR or move bases.

ANGSD has three BAQ options. BAQ = 0 turns it off, but the [documentation](https://www.popgen.dk/angsd/index.php/Input#Arguments) is [unclear](
https://github.com/ANGSD/angsd/issues/97#issuecomment-327461936) on what BAQ = 1 and BAQ = 2 actually do. The confusion probably comes from `samtools` using 'partial' to refer to heuristics for selecting *which* reads to realign and using 'extended' in older versions to refer to *how* the qualities were adjusted. 

In any event, the source code (annotations and any errors below are mine) shows that:
* -baq 1 ('partial' according to ANGSD but maybe 'simple' or 'local' are less misleading?): qualities for bases within each contiguous match/mismatch block (CIGAR M/=/X) that are flagged by the DP as an insertion, deletion, or shifted alignment are zeroed; the others are capped by the DP’s posterior Phred score.
* -baq 2 (extended): does the same zero/cap pass, then for each match block computes left-and-right running maxima over those provisional scores and takes the minimum at each position. This smooths low-confidence spikes into their neighbors, boosting sensitivity to SNPs near indels.

BAQ is expected to reduce false positive SNPs around indels ([Li 2011](https://doi.org/10.1093/bioinformatics/btr076](https://doi.org/10.1093/bioinformatics/btr076))). In practice, different BAQ models may lead to markedly different site frequency spectrum, but whether this is also a more accurate site frequency spectrum may be unclear and potentially dependent on factors such as sequencing depth (see discussions [here](https://www.biostars.org/p/9466154/), [here](https://github.com/ANGSD/angsd/issues/97) and [here](https://github.com/samtools/bcftools/pull/1363#issuecomment-80204203). 

```C
// Updated annotation including DP specifics
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

## Coefficient for capping mapping qualities (`-c`)

High mapping qualities do not necessarily indicate good alignments.  A read may align uniquely and earn MQ = 60 even if it carries dozens of high-quality mismatches (e.g. 20 × Q ≥ 20 → ΣQ ≥ 400) or large soft-clips (e.g. 50 × Q≈20 → ΣQ≈1,000).

The `-C` parameter is used to downgrade mapping qualities in cases like these. When invoked, the capped MQ function (`sam_cap_mapq` in the source code) first computes an error penalty *T* from the sum of high-quality mismatches (SubQ) and clipped-base penalties (ClipQ).  

Then, if *T* ≤ *C*, the new MQ is $= C×sqrt((C−T)/C)$ up to a maximum of *C*  and MQ = 0 if *T* > *C*. 

This down-weights moderately noisy alignments and eliminate the worst ones. Choosing *C* below the read aligner’s maximum MQ (e.g. < 60) will increasingly penalize true heterozygote reads (which carry one genuine mismatch), whereas higher values of  *C*  filter only the most egregious misalignments. 

## Usage

Provide example commands demonstrating a typical invocation:

```bash
bash <script1>.sh arg1 arg2
Rscript <script2>.R input_file output_file
```

---

## Inputs & Outputs

* **Inputs**:

  * `\<path/to/input1\>`: Description of the expected input file or directory.
  * `\<path/to/input2\>`: ...
* **Outputs**:

  * `\<path/to/output1\>`: Description of the generated output.
  * `\<path/to/output2\>`: ...

---

## Dependencies

List required modules, software, or packages:

* Bash (with `set -euo pipefail` recommended)
* [samtools](https://www.htslib.org/) >= 1.9
* R (packages: tidyverse, data.table, ...)
* ...

---

## Notes & Gotchas

* Any special instructions, known issues, or tips.
* For example: ensure you run this after the reduced reference is built.
* ...
