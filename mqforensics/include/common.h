#ifndef MQF_COMMON_H
#define MQF_COMMON_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>

#define BQ_MIN 13
#define BQ_CAP 33

// ---- Data types ----
typedef struct { char *chr; int start, end; } BedIv;

// small dynamic vector
typedef struct { double *v; int n, cap; } vecd;

// Per-site (analyze) accumulators
typedef struct {
    int depth;
    int mismatches;
    long long ins_len_sum, del_len_sum, clip_bases_sum;

    // direct-mode vectors (per-site)
    vecd mq, capmq, effmq, subq, clipq;

    // suffstats accumulators (used when --emit-suffstats)
    long long n_mq, n_cap, n_eff, n_subQ, n_clipQ, n_clipfrac;
    long double sum_mq, sum_cap, sum_eff, sum_subQ, sum_clipQ, sum_clipfrac;
    long double sumsq_mq, sumsq_cap, sumsq_eff, sumsq_subQ, sumsq_clipQ, sumsq_clipfrac;

    long long n_capped;             // count of reads where raw > cap
    long double sum_delta;          // S(raw - eff)
    long double sumsq_delta;        // S(raw - eff)^2

    // flanking context suffstats
    long long    n_flank_cov;
    long double  sum_flank_cov, sumsq_flank_cov;
    long long    n_flank_cf;
    long double  sum_flank_cf,  sumsq_flank_cf;

    // hist bins (per-BAM site; filled only when --emit-hist)
    int hist_mq[7];        // MQ [0..60] in bins of 10
    int hist_eff[7];       // effMQ [0..60] in bins of 10
    int hist_clipfrac[10]; // clip fraction [0..1] in bins of 0.1
} Site;

// Per-read QC (analyze)
typedef struct {
    int mapq_raw;
    int clipped_bases;   // soft+hard
    int M_q13;           // aligned bases with BQ>=13
    int X_ge13;          // mismatches with BQ>=13
    int SubQ;            // sum of min(bq,33) over mismatches
    int ClipQ;           // soft-clip qsum + 13*hardclip_len
    int cap_mq;          // raw cap in [0..C] or -1
    int cap_mq60;        // cap rescaled to [0..60]
    int eff_mq60;        // min(mapq_raw, cap_mq60)
    long long ins_len, del_len;
} ReadQC;

// Aggregator for summarize mode (pooled across BAMs)
typedef struct {
    char *chrom; int pos;
    long long depth;
    long long mismatch_bases, ins_len_sum, del_len_sum, clip_bases_sum;

    // moments: n, sum, sumsq
    long long n_mq, n_cap, n_eff, n_subQ, n_clipQ, n_clipfrac;
    long double sum_mq, sum_cap, sum_eff, sum_subQ, sum_clipQ, sum_clipfrac;
    long double sumsq_mq, sumsq_cap, sumsq_eff, sumsq_subQ, sumsq_clipQ, sumsq_clipfrac;

    // capping load
    long long n_capped;
    long double sum_delta, sumsq_delta; // raw - eff

    // NEW flanking suffstats
    long long    n_flank_cov;
    long double  sum_flank_cov, sumsq_flank_cov;
    long long    n_flank_cf;
    long double  sum_flank_cf,  sumsq_flank_cf;

    // pooled hists
    long long hist_mq[7], hist_eff[7], hist_clipfrac[10];
    bool have_hist;
} Agg;

// ---- Prototypes ----
// stats.c
void vpush(vecd *x, double a);
double dmean(const vecd *x);
double dmedian(vecd *x);
double dsd(const vecd *x);
double median_from_hist_uniform_bins(const long long *h, int nbins, double bin0_left, double bin_width);

// hist.c
void hist_accum_mq(int *hist7, int mq0_60);
void hist_accum_eff(int *hist7, int eff0_60);
void hist_accum_clipfrac(int *hist10, double frac0_1);
void hist_ks_w1_js(const long long *hA, const long long *hB, int nbins, double binw, double *ks, double *w1, double *js);

// nt16 ? uppercase base (A/C/G/T/N). Returns 'N' for ambiguous/15. 
static inline char nt16_to_base_uc(int nt16){
    static const char map[16] = {
        'N','A','C','N','G','N','N','N','T','N','N','N','N','N','N','N'
    };
    return map[(nt16 & 0xF)];
}

// bam_utils.c
int read_bed(const char *path, BedIv **out);
void free_bed(BedIv *a, int n);
void compute_read_qc(const bam1_t *b, int C, ReadQC *out);   // requires <htslib/sam.h>
long long apply_read_to_interval(const bam1_t *b, int iv_tid, int iv_start, int iv_end, const ReadQC *rq, Site *sites);
void add_local_mismatches(const bam1_t *b, int iv_tid, int iv_start, int iv_end, Site *sites);

// analyze.c
int analyze_main(int argc, char **argv);

// summarize.c
int summarize_main(int argc, char **argv);
void print_agg_header(FILE *fo, bool with_hist_metrics);
void finalize_and_print_agg(FILE *fo, Agg *a, bool with_hist_metrics);

#endif
