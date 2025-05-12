## Overview

Briefly describe the purpose of this directory (one or two sentences).

---

## Contents

* **`GATK_indeltargetcreator.sh`**: Short description of what this script does.
* **`GATK_indelrealigner.sh`**: Short description of what this analysis or step does.
---

## Important

* GATK v. >= 4 no longer includes the IndelTargetCreator and IndelRealigner walkers. 
* GATK v. 3 requires Java 8
* GATK requires a [sequence dictionary](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref) for indel realignment.
---

## Indel realignment

By default (O=6, E=1, B=4), BWA-MEM favors mismatches to single base-pair indels, and by convention, reports ungapped alignments over 2 bp indels. These unresolved small indels can result in false positive SNPs, and from previous experience with spruce, realignment reduces the number of SNPs with Ho >> 0.5 and high MAF but no observed heterozygotes. 

We performed indel realignment using GATK v. 3.7, which is implemented in two steps. First, [`RealignerTargetCreator`](https://github.com/broadgsa/gatk-protected/blob/master/public/gatk-tools-public/src/main/java/org/broadinstitute/gatk/tools/walkers/indels/RealignerTargetCreator.java#L171) creates a list of positions that could potentially benefit from realignment the CIGAR strings for each sample. Then, [`IndelRealigner`](https://github.com/broadgsa/gatk-protected/blob/master/public/gatk-tools-public/src/main/java/org/broadinstitute/gatk/tools/walkers/indels/IndelRealigner.java) performs the realignment at the identified target intervals.

## **`GATK_indeltargetcreator.sh`** usage

As input, we used the [intersected bams](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref) for each sample and limited the analysis to the regions in [`picea_newref_target_regions.bed`](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref) via the `-L` flag. 

```bash
#!/bin/bash
while IFS= read -r sample; do
    sbatch "$SCRIPTS/GATK_indeltargetcreator.sh" "$sample"
done < sample.list
```

## **`GATK_indelrealigner.sh`** usage

After `GATK_indeltargetcreator.sh` finishes, the target realignment intervals are saved to `${SPRUCE_PROJECT}/bams/realigned/realign_intervals` and supplied to `GATK_indelrealigner.sh`, which performs the realignment and outputs a new indexed BAM. We lowered `-entropy 0.05 [0.15]` and `-LOD 3 [5]` from their defaults to attempt realignment at more intervals.

```bash
#!/bin/bash
while IFS= read -r sample; do
    sbatch "$SCRIPTS/GATK_indelrealigner.sh" "$sample"
done < sample.list
```
---

## Inputs & Outputs

* **Inputs**:
  
  * `sample.list`: bam codes for samples that survived [intial QC](https://github.com/lxsllvn/spruceGBS/tree/main/03_initial_qc)
  * `${SPRUCE_PROJECT}/ref/picea_newref.fa`: the [reduced reference genome](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref)
  * `${SPRUCE_PROJECT}/ref/picea_newref_target_regions.bed`: the [target regions for analysis](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref) 
  * `${SPRUCE_PROJECT}/bams/intersected`: path to [bams intersected by `picea_newref_target_regions.bed`](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref) 

* **Outputs**:

  * `${SPRUCE_PROJECT}/bams/realigned/realign_intervals`: intervals identified for realignment by `GATK_indeltargetcreator.sh`; this is also a required input for `GATK_indelrealigner.sh`
  * `${SPRUCE_PROJECT}/bams/realigned`: realigned bams
---

## Dependencies

* [Java 8](https://github.com/adoptium/temurin8-binaries/releases); tested with jdk8u392-b08
* [GATK v. 3.7-0-gcfedb67](https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk;tab=objects?pli=1&prefix=&forceOnObjectsSortingFiltering=false). 
* [samtools](https://www.htslib.org/) v. 1.19.2
---


