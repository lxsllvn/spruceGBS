#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include "common.h"

static inline unsigned long long site_base_count(const Site *s, int slot){
    if (s->base_counts_hi){
        return s->base_counts_hi[slot];
    } else {
        const uint16_t *low = &s->nA_fwd;
        return low[slot];
    }
}

static void usage_an(void){
    fprintf(stderr,
      "Usage: mq_forensics analyze -b BAM -r BED -C INT -d MIN "
      "-o per_site.tsv -O per_interval.tsv "
      "[-f ref.fa] [--flank N] [--ref-only] "
      "[--emit-suffstats] [--emit-hist]\n");
}

static void write_analyze_headers(FILE *fo_site, int suffstats, int emit_hist, int flankN){
    if (!suffstats){
        fprintf(fo_site,
            "chrom\tpos\tdepth\tmismatch_bases\tins_len_sum\tdel_len_sum\tclip_bases_sum\t"
            "mq_mean\tmq_median\tmq_sd\t"
            "capmq60_mean\tcapmq60_median\tcapmq60_sd\t"
            "effmq_mean\teffmq_median\teffmq_sd\t"
            "subQ_mean\tsubQ_median\tsubQ_sd\t"
            "clipQ_mean\tclipQ_median\tclipQ_sd");
        if (flankN > 0){
            fprintf(fo_site, "\tflank_cov_mean\tflank_clipfrac_mean");
        }
        fprintf(fo_site,
            "\tnA_fwd\tnC_fwd\tnG_fwd\tnT_fwd\tnN_fwd"
            "\tnA_rev\tnC_rev\tnG_rev\tnT_rev\tnN_rev"
            "\tdepth_fwd\tdepth_rev"
            "\tentropy_pooled\talph_eff_pooled"
            "\tentropy_fwd\talph_eff_fwd"
            "\tentropy_rev\talph_eff_rev"
            "\tgc_frac_pooled\tgc_frac_fwd\tgc_frac_rev"
            "\tstrand_bias_z\n");
    } else {
        fprintf(fo_site,
            "chrom\tpos\tdepth\tmismatch_bases\tins_len_sum\tdel_len_sum\tclip_bases_sum\t"
            "n_mq\tsum_mq\tsumsq_mq\t"
            "n_cap\tsum_cap\tsumsq_cap\t"
            "n_eff\tsum_eff\tsumsq_eff\t"
            "n_subQ\tsum_subQ\tsumsq_subQ\t"
            "n_clipQ\tsum_clipQ\tsumsq_clipQ\t"
            "n_clipfrac\tsum_clipfrac\tsumsq_clipfrac\t"
            "n_capped\tsum_delta\tsumsq_delta\t"
            "n_flank_cov\tsum_flank_cov\tsumsq_flank_cov\t"
            "n_flank_cf\tsum_flank_cf\tsumsq_flank_cf\t"
            "nA_fwd\tnC_fwd\tnG_fwd\tnT_fwd\tnN_fwd\t"
            "nA_rev\tnC_rev\tnG_rev\tnT_rev\tnN_rev");
        if (emit_hist){
            for(int b=0;b<7;b++)  fprintf(fo_site, "\thist_mq_b%d", b);
            for(int b=0;b<7;b++)  fprintf(fo_site, "\thist_eff_b%d", b);
            for(int b=0;b<10;b++) fprintf(fo_site, "\thist_clipfrac_b%d", b);
        }
        fputc('\n', fo_site);
    }
}

int analyze_main(int argc, char **argv){
    const char *bam_path=NULL, *bed_path=NULL, *out_site=NULL, *out_iv=NULL;
    int C=50, min_depth=0, emit_suff=0, emit_hist=0;

    /* new flank/ref options */
    const char *fasta_path=NULL;
    int flankN = 0;
    int ref_only = 0;

    for (int i=1;i<argc;i++){
        if (!strcmp(argv[i],"-b") && i+1<argc) bam_path=argv[++i];
        else if (!strcmp(argv[i],"-r") && i+1<argc) bed_path=argv[++i];
        else if (!strcmp(argv[i],"-C") && i+1<argc) C=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-d") && i+1<argc) min_depth=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-o") && i+1<argc) out_site=argv[++i];
        else if (!strcmp(argv[i],"-O") && i+1<argc) out_iv=argv[++i];
        else if (!strcmp(argv[i],"--emit-suffstats")) emit_suff=1;
        else if (!strcmp(argv[i],"--emit-hist")) emit_hist=1;
        else if ( (!strcmp(argv[i],"-f") || !strcmp(argv[i],"--fasta")) && i+1<argc) fasta_path = argv[++i];
        else if (!strcmp(argv[i],"--flank") && i+1<argc) flankN = atoi(argv[++i]);
        else if (!strcmp(argv[i],"--ref-only")) ref_only = 1;
        else { usage_an(); return 2; }
    }
    if (!bam_path || !bed_path || !out_site || !out_iv){ usage_an(); return 2; }
    if ((flankN > 0 || ref_only) && !fasta_path){
        fprintf(stderr, "ERROR: --flank/--ref-only require -f/--fasta ref.fa\n");
        return 2;
    }

    samFile *fp = sam_open(bam_path,"r");
    if(!fp){ fprintf(stderr,"ERROR: open BAM %s\n", bam_path); return 1; }
    bam_hdr_t *hdr = sam_hdr_read(fp);
    if(!hdr){ fprintf(stderr,"ERROR: read header\n"); sam_close(fp); return 1; }
    hts_idx_t *idx = sam_index_load(fp, bam_path);
    if(!idx){ fprintf(stderr,"ERROR: BAM index not found for %s\n", bam_path); return 1; }

    BedIv *bed=NULL; int nbed = read_bed(bed_path, &bed);

    FILE *fo_site = fopen(out_site,"w");
    FILE *fo_iv   = fopen(out_iv,"w");
    if(!fo_site||!fo_iv){ fprintf(stderr,"ERROR: open outputs\n"); return 1; }

    write_analyze_headers(fo_site, emit_suff, emit_hist, flankN);
    fprintf(fo_iv, "chrom\tstart\tend\taligned_bp_in_interval\n");

    bam1_t *b=bam_init1();
    faidx_t *fai = NULL;
    if (fasta_path){
        fai = fai_load(fasta_path);
        if (!fai){ fprintf(stderr,"ERROR: cannot load FASTA %s\n", fasta_path); return 1; }
    }

    for (int iv=0; iv<nbed; iv++){
        int tid = sam_hdr_name2tid(hdr, bed[iv].chr);
        if (tid < 0){
            fprintf(stderr,"WARN: contig %s not in BAM header; skipping %s:%d-%d\n",
                    bed[iv].chr, bed[iv].chr, bed[iv].start, bed[iv].end);
            continue;
        }
        int L = bed[iv].end - bed[iv].start; if (L<=0) continue;
        Site *sites = (Site*)calloc(L, sizeof(Site)); if(!sites){ perror("calloc"); exit(1); }

        /* flanking accumulators (within the BED interval only) */
        long long *match_count = NULL, *cf_cnt = NULL;
        long double *cf_sum = NULL;
        char *refseq = NULL;

        if ((flankN > 0 || ref_only) && fai){
            match_count = (long long*)calloc(L, sizeof(long long));
            cf_cnt      = (long long*)calloc(L, sizeof(long long));
            cf_sum      = (long double*)calloc(L, sizeof(long double));
            if (!match_count || !cf_cnt || !cf_sum){ perror("calloc"); exit(1); }

            int seqlen = 0;
            refseq = faidx_fetch_seq(fai, bed[iv].chr, bed[iv].start, bed[iv].end-1, &seqlen);
            if (!refseq || seqlen != L){
                fprintf(stderr,"ERROR: faidx_fetch_seq(%s:%d-%d) failed (got %d)\n",
                        bed[iv].chr, bed[iv].start, bed[iv].end-1, seqlen);
                return 1;
            }
            /* uppercase ref */
            for (int k=0;k<L;k++){
                char c = refseq[k];
                if (c>='a' && c<='z') refseq[k]=c-32;
            }
        }

        hts_itr_t *itr = sam_itr_queryi(idx, tid, bed[iv].start, bed[iv].end);
        if (!itr){
            fprintf(stderr,"WARN: bad iterator for %s:%d-%d\n", bed[iv].chr, bed[iv].start, bed[iv].end);
            if (refseq) free(refseq);
            free(match_count); free(cf_cnt); free(cf_sum);
            free(sites); continue;
        }
        long long aligned_bp_sum = 0;

        while (sam_itr_next(fp, itr, b) >= 0){
            if (b->core.flag & (BAM_FUNMAP|BAM_FSECONDARY|BAM_FSUPPLEMENTARY)) continue;

            ReadQC rq; compute_read_qc(b, C, &rq);
            aligned_bp_sum += apply_read_to_interval(
                b, tid, bed[iv].start, bed[iv].end, &rq, sites
            );
            add_local_mismatches(b, tid, bed[iv].start, bed[iv].end, sites);

            /* If we’re doing suffstats/hist, attribute once per site (WR) */
            if (emit_suff || emit_hist){
                const uint32_t *cig=bam_get_cigar(b); int nc=b->core.n_cigar;
                int64_t ref = b->core.pos;
                for(int k=0;k<nc;k++){
                    int op=bam_cigar_op(cig[k]), ln=bam_cigar_oplen(cig[k]);
                    if (op==BAM_CMATCH||op==BAM_CEQUAL||op==BAM_CDIFF){
                        for(int i=0;i<ln;i++){
                            int64_t rp=ref+i;
                            if (rp>=bed[iv].start && rp<bed[iv].end){
                                int idxs = (int)(rp - bed[iv].start);
                                Site *s = &sites[idxs];

                                if (emit_suff){
                                    int raw=rq.mapq_raw, cap=rq.cap_mq60, eff=rq.eff_mq60;
                                    double d = (double)raw - (double)eff;
                                    s->n_mq++;    s->sum_mq += raw; s->sumsq_mq += (long double)raw*raw;
                                    s->n_cap++;   s->sum_cap += cap; s->sumsq_cap += (long double)cap*cap;
                                    s->n_eff++;   s->sum_eff += eff; s->sumsq_eff += (long double)eff*eff;
                                    s->n_subQ++;  s->sum_subQ += rq.SubQ; s->sumsq_subQ += (long double)rq.SubQ*rq.SubQ;
                                    s->n_clipQ++; s->sum_clipQ += rq.ClipQ; s->sumsq_clipQ += (long double)rq.ClipQ*rq.ClipQ;

                                    double cf = (double)rq.clipped_bases / (double)b->core.l_qseq;
                                    if (cf < 0) cf = 0;
                                    if (cf > 1) cf = 1;
                                    s->n_clipfrac++; s->sum_clipfrac += cf; s->sumsq_clipfrac += (long double)cf*cf;

                                    if (raw > cap) s->n_capped++;
                                    s->sum_delta += d; s->sumsq_delta += d*d;
                                }
                                if (emit_hist){
                                    hist_accum_mq(s->hist_mq, rq.mapq_raw);
                                    hist_accum_eff(s->hist_eff, rq.eff_mq60);
                                    double cf = (double)rq.clipped_bases / (double)b->core.l_qseq;
                                    if (cf < 0) cf = 0;
                                    if (cf > 1) cf = 1;
                                    hist_accum_clipfrac(s->hist_clipfrac, cf);
                                }
                            }
                        }
                        ref += ln;
                    } else if (op==BAM_CDEL || op==BAM_CREF_SKIP){
                        ref += ln;
                    }
                }
            }

            /* Flank accumulators (within BED bounds only) */
            if ((flankN > 0 || ref_only) && refseq){
                const uint32_t *cig=bam_get_cigar(b); int nc=b->core.n_cigar;
                int64_t rref = b->core.pos; int rpos=0;
                uint8_t *seq = bam_get_seq(b);
                /* clip fraction once per read (WR attribute) */
                double cf_read = (double)rq.clipped_bases / (double)b->core.l_qseq;
                for (int k=0;k<nc;k++){
                    int op=bam_cigar_op(cig[k]), ln=bam_cigar_oplen(cig[k]);
                    if (op==BAM_CMATCH||op==BAM_CEQUAL||op==BAM_CDIFF){
                        for (int i=0;i<ln;i++){
                            int64_t rp = rref + i;
                            if (rp >= bed[iv].start && rp < bed[iv].end){
                                int idxs = (int)(rp - bed[iv].start);
                                int nt = bam_seqi(seq, rpos+i) & 0xF;
                                char rb = nt16_to_base_uc(nt);
                                char rf = refseq[idxs];
                                int ok = 1;
                                if (ref_only){
                                    ok = (rb!='N' && rf!='N' && rb==rf);
                                }
                                if (ok){
                                    match_count[idxs] += 1;
                                    cf_sum[idxs] += cf_read;
                                    cf_cnt[idxs] += 1;
                                }
                            }
                        }
                        rref += ln; rpos += ln;
                    } else if (op==BAM_CINS){
                        rpos += ln;
                    } else if (op==BAM_CDEL || op==BAM_CREF_SKIP){
                        rref += ln;
                    } else if (op==BAM_CSOFT_CLIP){
                        rpos += ln;
                    }
                }
            }
        }
        hts_itr_destroy(itr);

        /* prefix sums for flank window means (for both direct and suffstats) */
        long long *pref_mc = NULL, *pref_cf_cnt = NULL;
        long double *pref_cf_sum = NULL;
        if (flankN > 0 && match_count){
            pref_mc     = (long long*)malloc((L+1)*sizeof(long long));
            pref_cf_cnt = (long long*)malloc((L+1)*sizeof(long long));
            pref_cf_sum = (long double*)malloc((L+1)*sizeof(long double));
            if (!pref_mc || !pref_cf_cnt || !pref_cf_sum){ perror("malloc"); exit(1); }
            pref_mc[0]=0; pref_cf_cnt[0]=0; pref_cf_sum[0]=0;
            for (int i=0;i<L;i++){
                pref_mc[i+1]     = pref_mc[i]     + match_count[i];
                pref_cf_cnt[i+1] = pref_cf_cnt[i] + cf_cnt[i];
                pref_cf_sum[i+1] = pref_cf_sum[i] + cf_sum[i];
            }
        }

        /* in suffstats mode, fold flank window means into per-site n/sum/sumsq */
        if (emit_suff && flankN > 0 && match_count){
            for (int i=0;i<L;i++){
                int left  = i - flankN; if (left  < 0) left  = 0;
                int right = i + flankN; if (right >= L) right = L - 1;
                int win_len = right - left + 1;

                long long mc_win       = pref_mc[right+1]     - pref_mc[left];
                long long cf_cnt_win   = pref_cf_cnt[right+1] - pref_cf_cnt[left];
                long double cf_sum_win = pref_cf_sum[right+1] - pref_cf_sum[left];

                double cov_mean = (double)mc_win / (double)win_len;
                double cf_mean  = (cf_cnt_win > 0) ? (double)(cf_sum_win / (long double)cf_cnt_win) : NAN;

                Site *s = &sites[i];
                /* coverage window mean always defined */
                s->n_flank_cov++;
                s->sum_flank_cov   += cov_mean;
                s->sumsq_flank_cov += (long double)cov_mean * cov_mean;

                if (!isnan(cf_mean)){
                    s->n_flank_cf++;
                    s->sum_flank_cf   += cf_mean;
                    s->sumsq_flank_cf += (long double)cf_mean * cf_mean;
                }
            }
        }

        /* write rows */
        for (int i = 0; i < L; i++) {
            int pos1 = bed[iv].start + i + 1;
            Site *s  = &sites[i];

            if (!emit_suff) {
                if (s->depth < min_depth) {
                    fprintf(fo_site, "%s\t%d", bed[iv].chr, pos1);
                    int na_fields = 20 + (flankN>0?2:0) + 22;
                    for (int c=0;c<na_fields;c++) fprintf(fo_site, "\tNA");
                    fputc('\n', fo_site);
                    if (s->base_counts_hi) {
                        free(s->base_counts_hi);
                        s->base_counts_hi = NULL;
                    }
                } else {
                    double mq_mu  = dmean(&s->mq);
                    double mq_med = dmedian(&s->mq);
                    double mq_sd  = dsd(&s->mq);

                    double cap_mu  = dmean(&s->capmq);
                    double cap_med = dmedian(&s->capmq);
                    double cap_sd  = dsd(&s->capmq);

                    double eff_mu  = dmean(&s->effmq);
                    double eff_med = dmedian(&s->effmq);
                    double eff_sd  = dsd(&s->effmq);

                    double sq_mu   = dmean(&s->subq);
                    double sq_med  = dmedian(&s->subq);
                    double sq_sd   = dsd(&s->subq);

                    double cq_mu   = dmean(&s->clipq);
                    double cq_med  = dmedian(&s->clipq);
                    double cq_sd   = dsd(&s->clipq);

                    double flank_cov_mean = NAN, flank_cf_mean = NAN;
                    if (flankN > 0 && pref_mc) {
                        int left  = i - flankN; if (left  < 0) left  = 0;
                        int right = i + flankN; if (right >= L) right = L - 1;
                        int win_len = right - left + 1;

                        long long     mc_win     = pref_mc[right + 1]     - pref_mc[left];
                        long long     cf_cnt_win = pref_cf_cnt[right + 1] - pref_cf_cnt[left];
                        long double   cf_sum_win = pref_cf_sum[right + 1] - pref_cf_sum[left];

                        flank_cov_mean = (double)mc_win / (double)win_len;
                        flank_cf_mean  = (cf_cnt_win > 0)
                                              ? (double)(cf_sum_win / (long double)cf_cnt_win)
                                              : NAN;
                    }

                    long long depth_fwd = s->nA_fwd + s->nC_fwd + s->nG_fwd + s->nT_fwd + s->nN_fwd;
                    long long depth_rev = s->nA_rev + s->nC_rev + s->nG_rev + s->nT_rev + s->nN_rev;

                    long long nA_tot = s->nA_fwd + s->nA_rev;
                    long long nC_tot = s->nC_fwd + s->nC_rev;
                    long long nG_tot = s->nG_fwd + s->nG_rev;
                    long long nT_tot = s->nT_fwd + s->nT_rev;

                    double ent_pooled = entropy_from_counts(nA_tot, nC_tot, nG_tot, nT_tot);
                    double ent_fwd    = entropy_from_counts(s->nA_fwd, s->nC_fwd, s->nG_fwd, s->nT_fwd);
                    double ent_rev    = entropy_from_counts(s->nA_rev, s->nC_rev, s->nG_rev, s->nT_rev);
                    double alph_p     = isnan(ent_pooled) ? NAN : pow(2.0, ent_pooled);
                    double alph_f     = isnan(ent_fwd)    ? NAN : pow(2.0, ent_fwd);
                    double alph_r     = isnan(ent_rev)    ? NAN : pow(2.0, ent_rev);
                    double gc_pooled  = gcfrac_from_counts(nA_tot, nC_tot, nG_tot, nT_tot);
                    double gc_fwd     = gcfrac_from_counts(s->nA_fwd, s->nC_fwd, s->nG_fwd, s->nT_fwd);
                    double gc_rev     = gcfrac_from_counts(s->nA_rev, s->nC_rev, s->nG_rev, s->nT_rev);
                    double sb_z       = strand_bias_z(depth_fwd, depth_rev);

                    fprintf(fo_site,
                            "%s\t%d\t%d\t%d\t%lld\t%lld\t%lld\t"
                            "%.3f\t%.3f\t%.3f\t"
                            "%.3f\t%.3f\t%.3f\t"
                            "%.3f\t%.3f\t%.3f\t"
                            "%.3f\t%.3f\t%.3f\t"
                            "%.3f\t%.3f\t%.3f",
                            bed[iv].chr, pos1,
                            s->depth, s->mismatches, s->ins_len_sum, s->del_len_sum, s->clip_bases_sum,
                            mq_mu,  mq_med,  mq_sd,
                            cap_mu, cap_med, cap_sd,
                            eff_mu, eff_med, eff_sd,
                            sq_mu,  sq_med,  sq_sd,
                            cq_mu,  cq_med,  cq_sd);
                    if (flankN > 0)
                        fprintf(fo_site, "\t%.6g\t%.6g", flank_cov_mean, flank_cf_mean);

                    fprintf(fo_site,
                            "\t%lld\t%lld\t%lld\t%lld\t%lld"
                            "\t%lld\t%lld\t%lld\t%lld\t%lld"
                            "\t%lld\t%lld"
                            "\t%.6g\t%.6g"
                            "\t%.6g\t%.6g"
                            "\t%.6g\t%.6g"
                            "\t%.6g\t%.6g\t%.6g"
                            "\t%.6g\n",
                            s->nA_fwd, s->nC_fwd, s->nG_fwd, s->nT_fwd, s->nN_fwd,
                            s->nA_rev, s->nC_rev, s->nG_rev, s->nT_rev, s->nN_rev,
                            depth_fwd, depth_rev,
                            ent_pooled, alph_p,
                            ent_fwd,   alph_f,
                            ent_rev,   alph_r,
                            gc_pooled, gc_fwd, gc_rev,
                            sb_z);
                }

                /* free vectors used in direct mode */
                free(s->mq.v);
                free(s->capmq.v);
                free(s->effmq.v);
                free(s->subq.v);
                free(s->clipq.v);

            } else {
                /* suffstats mode */
                if (s->depth < min_depth) {
                    fprintf(fo_site,
                        "%s\t%d\t%d\t%d\t%lld\t%lld\t%lld\t"
                        "0\t0\t0\t"  /* n_mq,     sum_mq,     sumsq_mq */
                        "0\t0\t0\t"  /* n_cap,    sum_cap,    sumsq_cap */
                        "0\t0\t0\t"  /* n_eff,    sum_eff,    sumsq_eff */
                        "0\t0\t0\t"  /* n_subQ,   sum_subQ,   sumsq_subQ */
                        "0\t0\t0\t"  /* n_clipQ,  sum_clipQ,  sumsq_clipQ */
                        "0\t0\t0\t"  /* n_clipfrac,sum_clipfrac,sumsq_clipfrac */
                        "0\t0\t0\t"  /* n_capped, sum_delta,  sumsq_delta */
                        "0\t0\t0\t"  /* n_flank_cov, sum_flank_cov, sumsq_flank_cov */
                        "0\t0\t0",   /* n_flank_cf,  sum_flank_cf,  sumsq_flank_cf */
                        bed[iv].chr, pos1,
                        s->depth, s->mismatches, s->ins_len_sum, s->del_len_sum, s->clip_bases_sum);
                    for (int k=0;k<10;k++) fputs("\t0", fo_site);
                    if (emit_hist) {
                        for (int b=0; b<7;  b++) fprintf(fo_site, "\t0");
                        for (int b=0; b<7;  b++) fprintf(fo_site, "\t0");
                        for (int b=0; b<10; b++) fprintf(fo_site, "\t0");
                    }
                    fputc('\n', fo_site);
                    if (s->base_counts_hi) {
                        free(s->base_counts_hi);
                        s->base_counts_hi = NULL;
                    }
                } else {
                    unsigned long long counts[10];
                    for (int cc=0; cc<10; cc++) counts[cc] = site_base_count(s, cc);
                    fprintf(fo_site,
                        "%s\t%d\t%d\t%d\t%lld\t%lld\t%lld\t"
                        "%lld\t%.10Lf\t%.10Lf\t"  /* mq */
                        "%lld\t%.10Lf\t%.10Lf\t"  /* cap */
                        "%lld\t%.10Lf\t%.10Lf\t"  /* eff */
                        "%lld\t%.10Lf\t%.10Lf\t"  /* subQ */
                        "%lld\t%.10Lf\t%.10Lf\t"  /* clipQ */
                        "%lld\t%.10Lf\t%.10Lf\t"  /* clipfrac */
                        "%lld\t%.10Lf\t%.10Lf\t"  /* capped/delta */
                        "%lld\t%.10Lf\t%.10Lf\t"  /* flank_cov */
                        "%lld\t%.10Lf\t%.10Lf\t"  /* flank_cf  */
                        "%lld\t%lld\t%lld\t%lld\t%lld\t"
                        "%lld\t%lld\t%lld\t%lld\t%lld",
                        bed[iv].chr, pos1,
                        s->depth, s->mismatches, s->ins_len_sum, s->del_len_sum, s->clip_bases_sum,
                        s->n_mq,         s->sum_mq,         s->sumsq_mq,
                        s->n_cap,        s->sum_cap,        s->sumsq_cap,
                        s->n_eff,        s->sum_eff,        s->sumsq_eff,
                        s->n_subQ,       s->sum_subQ,       s->sumsq_subQ,
                        s->n_clipQ,      s->sum_clipQ,      s->sumsq_clipQ,
                        s->n_clipfrac,   s->sum_clipfrac,   s->sumsq_clipfrac,
                        s->n_capped,     s->sum_delta,      s->sumsq_delta,
                        s->n_flank_cov,  s->sum_flank_cov,  s->sumsq_flank_cov,
                        s->n_flank_cf,   s->sum_flank_cf,   s->sumsq_flank_cf,
                        s->nA_fwd, s->nC_fwd, s->nG_fwd, s->nT_fwd, s->nN_fwd,
                        s->nA_rev, s->nC_rev, s->nG_rev, s->nT_rev, s->nN_rev);
                    if (emit_hist) {
                        for (int b=0; b<7;  b++) fprintf(fo_site, "\t%d", s->hist_mq[b]);
                        for (int b=0; b<7;  b++) fprintf(fo_site, "\t%d", s->hist_eff[b]);
                        for (int b=0; b<10; b++) fprintf(fo_site, "\t%d", s->hist_clipfrac[b]);
                    }
                    fputc('\n', fo_site);
                    if (s->base_counts_hi) {
                        free(s->base_counts_hi);
                        s->base_counts_hi = NULL;
                    }
                }
            }
        }

        fprintf(fo_iv, "%s\t%d\t%d\t%lld\n",
                bed[iv].chr, bed[iv].start, bed[iv].end, aligned_bp_sum);
        if (refseq) free(refseq);
        free(match_count); free(cf_cnt); free(cf_sum);
        free(pref_mc); free(pref_cf_cnt); free(pref_cf_sum);
        free(sites);
    }

    bam_destroy1(b);
    if (fai) fai_destroy(fai);
    fclose(fo_site); fclose(fo_iv);
    free_bed(bed, nbed);
    hts_idx_destroy(idx);
    bam_hdr_destroy(hdr);
    sam_close(fp);
    return 0;
}
