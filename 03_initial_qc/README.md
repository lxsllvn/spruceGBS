## Overview

Implements **Step 3: initial sample quality control** of the spruceGBS pipeline. This is a super-basic filter pass to remove samples with the lowest depth and/or breadth of coverage.  

---

## Contents

* **`bam_mapping_summary.sh`**: finds the number of mapped scaffolds and the total number of mapped reads per sample and returns a summary dataframe

---
## Sequencing depth and breadth per sample

While genotype likelihoods are useful for low-depth sequences, variation in depth and/or coverage among samples can bias estimates of diversity and structure.

In later steps, we systematically estimate the minimum depth/coverage required for analysis, but the lowest quality samples are not salvageable (I tried) and can be removed at this stage. We calculated the number of scaffolds with at least 1 mapped read and the total number of mapped reads for each sample using **bam_mapping_summary.sh**.

## **`bam_mapping_summary.sh`** usage

```bash
#!/bin/bash
sbatch "$SCRIPTS/bam_mapping_summary.sh" \
    "$SPRUCE_PROJECT/intersected_bam_mapping_summary" \ # outname
    "$SPRUCE_PROJECT/bams/intersected" # path to bams
```

## Visualization and removal

Samples had a median of 515,455 mapped reads (IQR = 762,618) to 7,983 (IQR = 2,724) scaffolds. Samples below the 10th percentile in either metrics (indicated by black, vertical lines) were removed (n = 208).

![sample_intial_qc](https://github.com/user-attachments/assets/24274798-021a-43ed-a290-872868966bf2)

Code to reproduce the figure and and list of failed samples is provided below. Note that `write.unix()` is a wrapper for read.table to enforce Unix line-endings on Windows, available here. 

```R
# Load table with bam_code, total_mapped_reads,
# and n_scaffolds_with_mapped_reads
depth_summary <- read.table("intersected_bam_mapping_summary.tsv",
                            header = TRUE)

# Find lowest 10%
cutoff_reads  <- quantile(log(depth_summary$n_reads),
                         probs = c(0.1))
cutoff_scaffs <- quantile(log(depth_summary$n_scaffolds),
                          probs = c(0.1))

# Plot of the distribution of read counts and scaffolds covered
library(ggplot2)

p <- ggplot(depth_summary) +
  geom_density(aes(x = log(n_reads)),
               fill = "dodgerblue", alpha = 0.4) +
  xlab("log(n reads)") +
  geom_vline(xintercept = cutoff_reads) +
  theme_minimal()

p1 <- ggplot(depth_summary) +
  geom_density(aes(x = log(n_scaffolds)),
               fill = "dodgerblue", alpha = 0.4) +
  xlab("log(n scaffolds)") +
  geom_vline(xintercept = cutoff_scaffs) +
  theme_minimal() +
  theme(axis.title.y = element_blank())

combined <- p + p1

ggsave("sample_initial_qc.pdf", plot = combined,
       width = 6.64, height = 3.07)

# Extract bam_codes of failed samples and create a file for them
fail <- subset(depth_summary,
               log(n_reads)     < cutoff_reads |
               log(n_scaffolds) < cutoff_scaffs)

df <- data.frame(
  bam_code    = fail$bam_code,
  n_reads     = fail$n_reads,
  n_scaffolds = fail$n_scaffolds,  # fixed typo
  flag        = "FAIL",
  reason      = "initial depth")

# Use write.table on Mac/Unix
write.table(df, "initial_qc_failed_samples.txt",
            row.names = FALSE, col.names = TRUE,
            quote = FALSE)
```

---

## Inputs & Outputs

* **Inputs**:
  * `$SPRUCE_PROJECT/bams/intersected`: BAMs after [intersection with target regions](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref).
    
* **Outputs**:
  * `$SPRUCE_PROJECT/intersected_bam_mapping_summary.tsv`: dataframe with the sample bam code (bam_code), number of mapped reads (n_reads), and  number of mapped scaffolds (n_scaffolds)
  * `intial_qc_failed_samples.txt`: dataframe of failed samples, with columns for bam_code, n_reads, n_scaffolds, flag, and the reason for the fail flag (reason)
  
---

## Dependencies

List required modules, software, or packages:

* [samtools](https://www.htslib.org/) v. 1.19.2
* R (packages: ggplot2)
* write.unix() for Windows

---
