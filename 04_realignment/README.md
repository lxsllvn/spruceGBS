## Overview

Implements **Step 4: indel realignment** of the spruceGBS pipeline using GATK 3.7.

---

## Contents

* **`GATK_indeltargetcreator.sh`**: Runs GATK RealignerTargetCreator to identify intervals likely to contain misaligned indels. Outputs .intervals files.
* **`GATK_indelrealigner.sh`**: Uses GATK IndelRealigner to locally realign reads across those intervals and produces new indexed bams.
  
---

## Important notes

* GATK v. ≥ 4 no longer includes IndelTargetCreator and IndelRealigner. These scripts require GATK v. 3 and were tested with GATK v. 3.7-0-gcfedb67
* GATK v. 3 requires Java 8 and was tested with AdoptOpenJDK 8u392
* A sequence dictionary for the reduced reference is required (see [02_reduced_ref](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref))

---

## Indel realignment

By default (O=6, E=1, B=4), BWA-MEM favors mismatches to single base-pair indels, and by convention, reports ungapped alignments over 2 bp indels. These unresolved small indels can result in false positive SNPs, and from previous experience with spruce, realignment reduces the number of SNPs with Ho >> 0.5 and high MAF but no observed heterozygotes. 

We performed indel realignment using GATK v. 3.7, which is implemented in two steps. First, [`RealignerTargetCreator`](https://github.com/broadgsa/gatk-protected/blob/master/public/gatk-tools-public/src/main/java/org/broadinstitute/gatk/tools/walkers/indels/RealignerTargetCreator.java#L171) scans the CIGAR strings of each bam and creates a list of `*.intervals` that could potentially benefit from realignment. Then, [`IndelRealigner`](https://github.com/broadgsa/gatk-protected/blob/master/public/gatk-tools-public/src/main/java/org/broadinstitute/gatk/tools/walkers/indels/IndelRealigner.java) realigns reads at the identified target intervals.

## **`GATK_indeltargetcreator.sh`** usage

As input, we used the [intersected bams](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref) for each sample and limited the analysis to the regions in [`picea_newref_target_regions.bed`](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref) via the `-L` flag. 

```bash
#!/bin/bash
while IFS= read -r sample; do
    sbatch "$SCRIPTS/GATK_indeltargetcreator.sh" "$sample"
done < sample.list
```

* **Inputs**:
    *  `sample.list`: list of bam codes of samples surviving [intial QC](https://github.com/lxsllvn/spruceGBS/tree/main/03_initial_qc)
     *  `${SPRUCE_PROJECT}/ref/picea_newref.fa`: the [reduced reference genome](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref)
     * `${SPRUCE_PROJECT}/ref/picea_newref_target_regions.bed`: the [target regions for analysis](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref) 
      * `${SPRUCE_PROJECT}/bams/intersected`: path to [bams intersected by `picea_newref_target_regions.bed`](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref) 

* **Outputs**:
    * `${SPRUCE_PROJECT}/bams/realigned/realign_intervals/*.intervals`: intervals for realignment 

## **`GATK_indelrealigner.sh`** usage

Once the target realignment intervals are identified, `GATK_indelrealigner.sh` performs the realignment and outputs a new indexed BAM. We lowered `-entropy 0.05 [0.15]` and `-LOD 3 [5]` from their defaults to attempt realignment at more intervals.

```bash
#!/bin/bash
while IFS= read -r sample; do
    sbatch "$SCRIPTS/GATK_indelrealigner.sh" "$sample"
done < sample.list
```

* **Inputs**:
    *  `sample.list`: list of bam codes of samples surviving [intial QC]
    *  `${SPRUCE_PROJECT}/ref/picea_newref.fa`: the [reduced reference genome](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref)
    * `${SPRUCE_PROJECT}/bams/intersected`: path to [intersected bams](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref) 
    * `${SPRUCE_PROJECT}/bams/realigned/realign_intervals/*.intervals`: intervals for realignment

* **Outputs**:
 * `${SPRUCE_PROJECT}/bams/realigned/*.realigned.bam`: realigned, indexed bams

---

## Dependencies

* [Java 8](https://github.com/adoptium/temurin8-binaries/releases); tested with jdk8u392-b08
* [GATK v. 3.7-0-gcfedb67](https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk;tab=objects?pli=1&prefix=&forceOnObjectsSortingFiltering=false). 
* [samtools](https://www.htslib.org/) v. 1.19.2

---


