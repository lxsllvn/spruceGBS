#define _GNU_SOURCE
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "common.h"

static void usage_sum(void){
    fprintf(stderr, "Usage: mq_forensics summarize -i suffstats.sorted.tsv -o pooled.tsv\n");
}

static char* xstrdup_local(const char* s){
    size_t n = strlen(s)+1;
    char *p = (char*)malloc(n);
    if(!p){ perror("malloc"); exit(1); }
    memcpy(p,s,n);
    return p;
}

void print_agg_header(FILE *fo, bool with_hist_metrics){
    fprintf(fo,
        "chrom\tpos\tdepth\t"
        "mismatch_rate\tins_rate\tdel_rate\tclip_rate\t"
        "mq_mean\tmq_sd\tcapmq60_mean\tcapmq60_sd\t"
        "effmq_mean\teffmq_sd\t"
        "subQ_mean\tsubQ_sd\t"
        "clipQ_mean\tclipQ_sd\t"
        "clipfrac_mean\tclipfrac_sd\t"
        "frac_capped\tmean_delta\tsd_delta\t"
        "flank_cov_mean\tflank_cov_sd\t"
        "flank_cf_mean\tflank_cf_sd\t"
        "entropy_pooled\talph_eff_pooled\tentropy_fwd\talph_eff_fwd\t"
        "entropy_rev\talph_eff_rev\t"
        "gc_frac_pooled\tgc_frac_fwd\tgc_frac_rev\t"
        "strand_bias_z");
    if (with_hist_metrics){
        fprintf(fo, "\tks_mq_eff\tw1_mq_eff\tjs_mq_eff\tmq_median_hist\teffmq_median_hist\tclipfrac_median_hist");
    }
    fputc('\n', fo);
}
void finalize_and_print_agg(FILE *fo, Agg *a, bool with_hist_metrics){
    if (!a->chrom) return;

    // ---- pooled means & SDs (population variance from sums/sumsq) ----
    double mq_mean    = (a->n_mq>0)       ? (double)(a->sum_mq / a->n_mq) : NAN;
    double mq_varp    = (a->n_mq>0)       ? (double)(a->sumsq_mq / a->n_mq - (a->sum_mq/a->n_mq)*(a->sum_mq/a->n_mq)) : NAN;
    if (mq_varp < 0) mq_varp = 0;
    double mq_sd      = (a->n_mq>0)       ? sqrt(mq_varp) : NAN;

    double cap_mean   = (a->n_cap>0)      ? (double)(a->sum_cap / a->n_cap) : NAN;
    double cap_varp   = (a->n_cap>0)      ? (double)(a->sumsq_cap / a->n_cap - (a->sum_cap/a->n_cap)*(a->sum_cap/a->n_cap)) : NAN;
    if (cap_varp < 0) cap_varp = 0;
    double cap_sd     = (a->n_cap>0)      ? sqrt(cap_varp) : NAN;

    double eff_mean   = (a->n_eff>0)      ? (double)(a->sum_eff / a->n_eff) : NAN;
    double eff_varp   = (a->n_eff>0)      ? (double)(a->sumsq_eff / a->n_eff - (a->sum_eff/a->n_eff)*(a->sum_eff/a->n_eff)) : NAN;
    if (eff_varp < 0) eff_varp = 0;
    double eff_sd     = (a->n_eff>0)      ? sqrt(eff_varp) : NAN;

    double subQ_mean  = (a->n_subQ>0)     ? (double)(a->sum_subQ / a->n_subQ) : NAN;
    double subQ_varp  = (a->n_subQ>0)     ? (double)(a->sumsq_subQ / a->n_subQ - (a->sum_subQ/a->n_subQ)*(a->sum_subQ/a->n_subQ)) : NAN;
    if (subQ_varp < 0) subQ_varp = 0;
    double subQ_sd    = (a->n_subQ>0)     ? sqrt(subQ_varp) : NAN;

    double clipQ_mean = (a->n_clipQ>0)    ? (double)(a->sum_clipQ / a->n_clipQ) : NAN;
    double clipQ_varp = (a->n_clipQ>0)    ? (double)(a->sumsq_clipQ / a->n_clipQ - (a->sum_clipQ/a->n_clipQ)*(a->sum_clipQ/a->n_clipQ)) : NAN;
    if (clipQ_varp < 0) clipQ_varp = 0;
    double clipQ_sd   = (a->n_clipQ>0)    ? sqrt(clipQ_varp) : NAN;

    double cf_mean    = (a->n_clipfrac>0) ? (double)(a->sum_clipfrac / a->n_clipfrac) : NAN;
    double cf_varp    = (a->n_clipfrac>0) ? (double)(a->sumsq_clipfrac / a->n_clipfrac - (a->sum_clipfrac/a->n_clipfrac)*(a->sum_clipfrac/a->n_clipfrac)) : NAN;
    if (cf_varp < 0) cf_varp = 0;
    double cf_sd      = (a->n_clipfrac>0) ? sqrt(cf_varp) : NAN;

    // capping deltas
    double frac_capped = (a->depth>0) ? (double)a->n_capped / (double)a->depth : NAN;
    double mean_delta  = (a->depth>0) ? (double)(a->sum_delta / (long double)a->depth) : NAN;
    double var_delta   = (a->depth>0) ? (double)(a->sumsq_delta / (long double)a->depth
                                 - (a->sum_delta/a->depth)*(a->sum_delta/a->depth)) : NAN;
    if (var_delta < 0) var_delta = 0;
    double sd_delta    = (a->depth>0) ? sqrt(var_delta) : NAN;

    // pooled global rates
    double mismatch_rate = (a->depth>0) ? (double)a->mismatch_bases / (double)a->depth : NAN;
    double ins_rate      = (a->depth>0) ? (double)a->ins_len_sum    / (double)a->depth : NAN;
    double del_rate      = (a->depth>0) ? (double)a->del_len_sum    / (double)a->depth : NAN;
    double clip_rate     = (a->depth>0) ? (double)a->clip_bases_sum / (double)a->depth : NAN;

    // divergences + histogram medians (if we have hist)
    double ks=NAN, w1=NAN, js=NAN;
    double mq_med_hist=NAN, eff_med_hist=NAN, cf_med_hist=NAN;
    if (with_hist_metrics && a->have_hist){
        hist_ks_w1_js(a->hist_mq, a->hist_eff, 7, 10.0, &ks, &w1, &js);
        mq_med_hist  = median_from_hist_uniform_bins(a->hist_mq,  7, 0.0, 10.0);
        eff_med_hist = median_from_hist_uniform_bins(a->hist_eff, 7, 0.0, 10.0);
        if (mq_med_hist  > 60.0) mq_med_hist  = 60.0;
        if (eff_med_hist > 60.0) eff_med_hist = 60.0;
        long long tot_cf=0; for(int b=0;b<10;b++) tot_cf+=a->hist_clipfrac[b];
        if (tot_cf > 0){
            cf_med_hist = median_from_hist_uniform_bins(a->hist_clipfrac, 10, 0.0, 0.1);
            if (cf_med_hist < 0.0) cf_med_hist = 0.0;
            if (cf_med_hist > 1.0) cf_med_hist = 1.0;
        }
    }

    /* --- pooled flanking stats (compute before printing) --- */
    double flank_cov_mean = (a->n_flank_cov > 0) ? (double)(a->sum_flank_cov / a->n_flank_cov) : NAN;
    double flank_cov_varp = (a->n_flank_cov > 0) ? (double)(a->sumsq_flank_cov / a->n_flank_cov
                                        - (a->sum_flank_cov/a->n_flank_cov)*(a->sum_flank_cov/a->n_flank_cov)) : NAN;
    if (flank_cov_varp < 0) flank_cov_varp = 0;
    double flank_cov_sd   = (a->n_flank_cov > 0) ? sqrt(flank_cov_varp) : NAN;

    double flank_cf_mean  = (a->n_flank_cf > 0) ? (double)(a->sum_flank_cf / a->n_flank_cf) : NAN;
    double flank_cf_varp  = (a->n_flank_cf > 0) ? (double)(a->sumsq_flank_cf / a->n_flank_cf
                                        - (a->sum_flank_cf/a->n_flank_cf)*(a->sum_flank_cf/a->n_flank_cf)) : NAN;
    if (flank_cf_varp < 0) flank_cf_varp = 0;
    double flank_cf_sd    = (a->n_flank_cf > 0) ? sqrt(flank_cf_varp) : NAN;

    long long depth_fwd = a->nA_fwd + a->nC_fwd + a->nG_fwd + a->nT_fwd + a->nN_fwd;
    long long depth_rev = a->nA_rev + a->nC_rev + a->nG_rev + a->nT_rev + a->nN_rev;
    long long nA_tot = a->nA_fwd + a->nA_rev;
    long long nC_tot = a->nC_fwd + a->nC_rev;
    long long nG_tot = a->nG_fwd + a->nG_rev;
    long long nT_tot = a->nT_fwd + a->nT_rev;
    double ent_pooled = entropy_from_counts(nA_tot, nC_tot, nG_tot, nT_tot);
    double ent_fwd    = entropy_from_counts(a->nA_fwd, a->nC_fwd, a->nG_fwd, a->nT_fwd);
    double ent_rev    = entropy_from_counts(a->nA_rev, a->nC_rev, a->nG_rev, a->nT_rev);
    double alph_p     = isnan(ent_pooled) ? NAN : pow(2.0, ent_pooled);
    double alph_f     = isnan(ent_fwd)    ? NAN : pow(2.0, ent_fwd);
    double alph_r     = isnan(ent_rev)    ? NAN : pow(2.0, ent_rev);
    double gc_pooled  = gcfrac_from_counts(nA_tot, nC_tot, nG_tot, nT_tot);
    double gc_fwd     = gcfrac_from_counts(a->nA_fwd, a->nC_fwd, a->nG_fwd, a->nT_fwd);
    double gc_rev     = gcfrac_from_counts(a->nA_rev, a->nC_rev, a->nG_rev, a->nT_rev);
    double sb_z       = strand_bias_z(depth_fwd, depth_rev);

    /* --- print exactly one row --- */
    fprintf(fo,
        "%s\t%d\t%lld\t"
        "%.6g\t%.6g\t%.6g\t%.6g\t"          /* mismatch/ins/del/clip rates */
        "%.6g\t%.6g\t%.6g\t%.6g\t"          /* mq mean/sd, cap mean/sd     */
        "%.6g\t%.6g\t%.6g\t%.6g\t"          /* eff mean/sd, subQ mean/sd   */
        "%.6g\t%.6g\t"                      /* clipQ mean/sd                */
        "%.6g\t%.6g\t"                      /* clipfrac mean/sd             */
        "%.6g\t%.6g\t%.6g\t"                /* frac_capped, mean_delta, sd_delta */
        "%.6g\t%.6g\t"                      /* flank_cov mean/sd            */
        "%.6g\t%.6g\t"                      /* flank_cf mean/sd             */
        "%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t"  /* entropy/alphabet sizes */
        "%.6g\t%.6g\t%.6g\t"                  /* GC fractions */
        "%.6g",                               /* strand bias z */
        a->chrom, a->pos, a->depth,
        mismatch_rate, ins_rate, del_rate, clip_rate,
        mq_mean, mq_sd, cap_mean, cap_sd,
        eff_mean, eff_sd, subQ_mean, subQ_sd,
        clipQ_mean, clipQ_sd,
        cf_mean, cf_sd,
        frac_capped, mean_delta, sd_delta,
        flank_cov_mean, flank_cov_sd,
        flank_cf_mean,  flank_cf_sd,
        ent_pooled, alph_p, ent_fwd, alph_f, ent_rev, alph_r,
        gc_pooled, gc_fwd, gc_rev,
        sb_z
    );

    if (with_hist_metrics){
        fprintf(fo, "\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g",
                ks, w1, js, mq_med_hist, eff_med_hist, cf_med_hist);
    }
    fputc('\n', fo);
}

static void zero_agg(Agg *a){
    memset(a, 0, sizeof(*a));
    a->chrom = NULL;
}

int summarize_main(int argc, char **argv){
    const char *in_path=NULL, *out_path=NULL;
    for (int i=1;i<argc;i++){
        if (!strcmp(argv[i],"-i") && i+1<argc) in_path=argv[++i];
        else if (!strcmp(argv[i],"-o") && i+1<argc) out_path=argv[++i];
        else { usage_sum(); return 2; }
    }
    if (!in_path || !out_path){ usage_sum(); return 2; }

    FILE *fi = strcmp(in_path,"-")==0? stdin : fopen(in_path,"r");
    if(!fi){ fprintf(stderr,"ERROR: open input %s\n", in_path); return 1; }
    FILE *fo = fopen(out_path,"w");
    if(!fo){ fprintf(stderr,"ERROR: open output %s\n", out_path); return 1; }

    char *line=NULL; size_t cap=0;
    ssize_t len;

    // Detect header & histogram presence (look for hist_mq_b0 token)
    int has_header=0; bool with_hist=false;
    {
        off_t pos = ftello(fi);
        char buf[8192];
        if (fgets(buf, sizeof(buf), fi)){
            if (!strncmp(buf,"chrom\tpos\t",11)) has_header=1;
            with_hist = (strstr(buf, "\thist_mq_b0") != NULL);
            fseeko(fi, pos, SEEK_SET);
        }
    }

    print_agg_header(fo, with_hist);
    if (has_header) getline(&line, &cap, fi); // consume header

    Agg cur; zero_agg(&cur);
    char *cur_chr=NULL; int cur_pos=-1;

    while ((len = getline(&line, &cap, fi)) != -1){
        if (len==0) continue;
        if (line[len-1]=='\n') line[len-1]=0;

        char *p = line, *tok;

        // chrom
        tok = strsep(&p, "\t"); if (!tok) continue;
        char *chrom = tok;
        // pos
        tok = strsep(&p, "\t"); if (!tok) continue;
        int pos = atoi(tok);

        // new group?
        if (!cur.chrom || strcmp(chrom, cur_chr) || pos != cur_pos){
            if (cur.chrom){
                finalize_and_print_agg(fo, &cur, with_hist);
                free(cur.chrom);
            }
            zero_agg(&cur);
            cur.chrom = xstrdup_local(chrom);
            cur_chr   = cur.chrom;
            cur_pos   = pos;
            cur.pos   = pos;
        }

        // depth, mismatch_bases, ins_len_sum, del_len_sum, clip_bases_sum
        tok = strsep(&p,"\t"); if(!tok) continue; long long depth = atoll(tok);
        tok = strsep(&p,"\t"); if(!tok) continue; long long mm    = atoll(tok);
        tok = strsep(&p,"\t"); if(!tok) continue; long long ins   = atoll(tok);
        tok = strsep(&p,"\t"); if(!tok) continue; long long del   = atoll(tok);
        tok = strsep(&p,"\t"); if(!tok) continue; long long clips = atoll(tok);

        cur.depth          += depth;
        cur.mismatch_bases += mm;
        cur.ins_len_sum    += ins;
        cur.del_len_sum    += del;
        cur.clip_bases_sum += clips;

        // six blocks: (n,sum,sumsq) x 6
        long long *n_ptr[6]      = { &cur.n_mq, &cur.n_cap, &cur.n_eff, &cur.n_subQ, &cur.n_clipQ, &cur.n_clipfrac };
        long double *s_ptr[6]    = { &cur.sum_mq, &cur.sum_cap, &cur.sum_eff, &cur.sum_subQ, &cur.sum_clipQ, &cur.sum_clipfrac };
        long double *ss_ptr[6]   = { &cur.sumsq_mq, &cur.sumsq_cap, &cur.sumsq_eff, &cur.sumsq_subQ, &cur.sumsq_clipQ, &cur.sumsq_clipfrac };

        for (int blk=0; blk<6; blk++){
            tok = strsep(&p,"\t"); if(!tok) goto parse_error; long long n = atoll(tok);
            tok = strsep(&p,"\t"); if(!tok) goto parse_error; long double s;  { char *end; s  = strtold(tok,&end); }
            tok = strsep(&p,"\t"); if(!tok) goto parse_error; long double ss; { char *end; ss = strtold(tok,&end); }
            *n_ptr[blk]  += n;
            *s_ptr[blk]  += s;
            *ss_ptr[blk] += ss;
        }

        // capping counters
        tok = strsep(&p,"\t"); if(!tok) goto parse_error; long long ncap = atoll(tok);
        tok = strsep(&p,"\t"); if(!tok) goto parse_error; long double sdel;  { char *end; sdel  = strtold(tok,&end); }
        tok = strsep(&p,"\t"); if(!tok) goto parse_error; long double ssdel; { char *end; ssdel = strtold(tok,&end); }
        cur.n_capped   += ncap;
        cur.sum_delta  += sdel;
        cur.sumsq_delta+= ssdel;

        // optional flanking suffstats (2 triplets). If absent, they stay zero.
        // n_flank_cov, sum_flank_cov, sumsq_flank_cov
        if (p && *p){
            tok = strsep(&p,"\t");
            if (tok){ long long n = atoll(tok);
                tok = strsep(&p,"\t"); if(!tok) goto parse_error; long double s;  { char *end; s  = strtold(tok,&end); }
                tok = strsep(&p,"\t"); if(!tok) goto parse_error; long double ss; { char *end; ss = strtold(tok,&end); }
                cur.n_flank_cov     += n;
                cur.sum_flank_cov   += s;
                cur.sumsq_flank_cov += ss;
            } else goto parse_error;
        }

        // n_flank_cf, sum_flank_cf, sumsq_flank_cf
        if (p && *p){
            tok = strsep(&p,"\t");
            if (tok){ long long n = atoll(tok);
                tok = strsep(&p,"\t"); if(!tok) goto parse_error; long double s;  { char *end; s  = strtold(tok,&end); }
                tok = strsep(&p,"\t"); if(!tok) goto parse_error; long double ss; { char *end; ss = strtold(tok,&end); }
                cur.n_flank_cf     += n;
                cur.sum_flank_cf   += s;
                cur.sumsq_flank_cf += ss;
            } else goto parse_error;
        }

        // per-strand base counts (10 numbers always present)
        long long *cnt_ptr[10] = { &cur.nA_fwd, &cur.nC_fwd, &cur.nG_fwd, &cur.nT_fwd, &cur.nN_fwd,
                                   &cur.nA_rev, &cur.nC_rev, &cur.nG_rev, &cur.nT_rev, &cur.nN_rev };
        for (int cc=0; cc<10; cc++){ tok = strsep(&p,"\t"); if(!tok) goto parse_error; *(cnt_ptr[cc]) += atoll(tok); }

        // optional hists
        if (with_hist){
            for (int b=0;b<7;b++){ tok=strsep(&p,"\t"); if(!tok) goto parse_error; cur.hist_mq[b]        += atoll(tok); }
            for (int b=0;b<7;b++){ tok=strsep(&p,"\t"); if(!tok) goto parse_error; cur.hist_eff[b]       += atoll(tok); }
            for (int b=0;b<10;b++){ tok=strsep(&p,"\t"); if(!tok) goto parse_error; cur.hist_clipfrac[b] += atoll(tok); }
            cur.have_hist = true;
        }
        continue;

parse_error:
        fprintf(stderr, "WARN: malformed line at %s:%d — skipping\n", chrom, pos);
    }

    if (cur.chrom){
        finalize_and_print_agg(fo, &cur, with_hist);
        free(cur.chrom);
    }

    free(line);
    fclose(fi); fclose(fo);
    return 0;
}
