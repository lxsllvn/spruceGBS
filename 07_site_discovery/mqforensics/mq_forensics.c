// mq_forensics.c
// Per-site (WR) and per-interval stats from a BAM over BED regions.
//   - Per site: depth, mismatches (local, from MD), INS/DEL length sums (local),
//               sum of clipped bases (WR), MAPQ mean/median/SD (raw),
//               capped MAPQ mean/median/SD (WR, with -C) rescaled to 0..60,
//               effective MAPQ = min(raw, capped) (WR, 0..60),
//               SubQ mean/median/SD (WR), ClipQ mean/median/SD (WR).
//   - Per interval: aligned_bp_in_interval = Σ overlap of aligned bases (M/= /X)
//                   with the interval (area under coverage).
//
// Build (with an htslib module loaded):
//   gcc -O3 -std=c11 mq_forensics.c -o mq_forensics $(pkg-config --cflags --libs htslib) -lm
//
// Usage:
//   ./mq_forensics -b sample.bam -r regions.bed -C 50 -d 0 \
//     -o sample.per_site.tsv -O sample.per_interval.tsv
//
// Notes:
//   - WR attribution: each read’s {MAPQ, capMQ60, effMQ60, SubQ, ClipQ, clipped_bases}
//     is added once to *every site* it overlaps in the region.
//   - Local counts (mismatches, INS/DEL starts) follow mpileup semantics.
//   - Requires MD/NM. If absent, run once: samtools calmd -bAr ref.fa in.bam > out.bam

#include <htslib/sam.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#define BQ_MIN 13
#define BQ_CAP 33

// portable strdup (some toolchains gate strdup behind POSIX macros)
static char* xstrdup(const char* s){
    size_t n = strlen(s) + 1;
    char *p = (char*)malloc(n);
    if (!p) { perror("malloc"); exit(1); }
    memcpy(p, s, n);
    return p;
}

// simple dynamic vector of doubles for WR per-site stats
typedef struct { double *v; int n, cap; } vecd;
static inline void vpush(vecd *x, double a){
    if(x->n==x->cap){
        x->cap = x->cap? x->cap*2:32;
        x->v = (double*)realloc(x->v, x->cap*sizeof(double));
        if(!x->v){ perror("realloc"); exit(1); }
    }
    x->v[x->n++]=a;
}
static int dcmp(const void*a,const void*b){
    double x=*(double*)a, y=*(double*)b;
    return (x<y)?-1:(x>y);
}
static inline double dmean(const vecd *x){
    if(!x->n) return 0.0;
    double s=0; for(int i=0;i<x->n;i++) s+=x->v[i]; return s/x->n;
}
static inline double dmedian(vecd *x){
    if(!x->n) return 0.0;
    qsort(x->v,x->n,sizeof(double),dcmp);
    if(x->n&1) return x->v[x->n/2];
    return 0.5*(x->v[x->n/2-1]+x->v[x->n/2]);
}
static inline double dsd(const vecd *x){
    if(x->n<2) return 0.0;
    double mu=dmean(x), s2=0;
    for(int i=0;i<x->n;i++){ double d=x->v[i]-mu; s2+=d*d; }
    return sqrt(s2/(x->n-1));
}

typedef struct {
    // per-site accumulators
    int depth;                 // includes reads with deletion at the site
    int mismatches;            // local mismatches at this ref pos (from MD)
    long long ins_len_sum;     // sum of insertion lengths starting here
    long long del_len_sum;     // sum of deletion lengths starting here
    long long clip_bases_sum;  // WR: soft+hard clipped bases of overlapping reads
    vecd mq;                   // WR: raw MAPQ per overlapping read (0..60)
    vecd capmq;                // WR: capped MAPQ rescaled to 0..60
    vecd effmq;                // WR: effective MQ = min(raw, capmq60) (0..60)
    vecd subq;                 // WR: SubQ per overlapping read
    vecd clipq;                // WR: ClipQ per overlapping read
} Site;

typedef struct { char *chr; int start, end; } BedIv;

// --- BAM helpers ---
static inline int left_clip_len(const bam1_t *b, int op){
    const uint32_t *cig=bam_get_cigar(b);
    if (b->core.n_cigar==0) return 0;
    return (bam_cigar_op(cig[0])==op) ? bam_cigar_oplen(cig[0]) : 0;
}
static inline int right_clip_len(const bam1_t *b, int op){
    const uint32_t *cig=bam_get_cigar(b); int n=b->core.n_cigar; if(!n) return 0;
    return (bam_cigar_op(cig[n-1])==op) ? bam_cigar_oplen(cig[n-1]) : 0;
}

// Mismatch record
typedef struct { int rp; int bq; } Mismatch;

// Extract mismatches via MD; returns array (caller frees) and count in *n_out
static Mismatch* mismatches_from_MD(const bam1_t *b, int *n_out){
    uint8_t *mdp = bam_aux_get(b,"MD");
    if (!mdp){ *n_out=0; return NULL; }
    const char *md = bam_aux2Z(mdp); if(!md){ *n_out=0; return NULL; }

    // Build aligned read positions: indices (in read) that align to reference (M/= /X)
    const uint32_t *cig = bam_get_cigar(b);
    int n_cigar = b->core.n_cigar;
    int rpos=0, cap=0;
    for(int k=0;k<n_cigar;k++){
        int op=bam_cigar_op(cig[k]), ln=bam_cigar_oplen(cig[k]);
        if (op==BAM_CMATCH||op==BAM_CEQUAL||op==BAM_CDIFF) cap+=ln;
        else if (op==BAM_CINS||op==BAM_CSOFT_CLIP) rpos+=ln;
    }
    int *aligned = cap? (int*)malloc(cap*sizeof(int)) : NULL;
    if (cap && !aligned){ perror("malloc"); exit(1); }
    int idx=0; rpos=0;
    for(int k=0;k<n_cigar;k++){
        int op=bam_cigar_op(cig[k]), ln=bam_cigar_oplen(cig[k]);
        if (op==BAM_CMATCH||op==BAM_CEQUAL||op==BAM_CDIFF){
            for(int i=0;i<ln;i++) aligned[idx++]=rpos++;
        } else if (op==BAM_CINS||op==BAM_CSOFT_CLIP){
            rpos+=ln;
        }
    }
    uint8_t *qual = bam_get_qual(b);
    int ai=0; const char *p=md;
    Mismatch *mm = cap? (Mismatch*)malloc(cap*sizeof(Mismatch)) : NULL; // upper bound
    if (cap && !mm){ perror("malloc"); exit(1); }
    int mmn=0;

    while (*p){
        if (isdigit((unsigned char)*p)){
            int v=0; while (isdigit((unsigned char)*p)){ v = v*10 + (*p-'0'); ++p; }
            ai += v;
        } else if (*p=='^'){
            ++p; while (*p && isalpha((unsigned char)*p)) ++p; // skip deletion
        } else if (isalpha((unsigned char)*p)){
            if (ai < cap){
                int rp = aligned[ai];
                int bq = qual ? qual[rp] : 0;
                mm[mmn].rp = rp; mm[mmn].bq = bq; mmn++;
                ai++;
            }
            ++p;
        } else ++p;
    }
    free(aligned);
    *n_out = mmn;
    return mm;
}

// Per-read metrics for WR attribution and capped MQ
typedef struct {
    int mapq_raw;        // BAM MAPQ (0..60)
    int clipped_bases;   // soft+hard
    int M_q13;           // # aligned bases with baseQ>=13 (M/= /X)
    int X_ge13;          // # mismatches with baseQ>=13
    int SubQ;            // sum of min(bq,33) over those mismatches
    int ClipQ;           // soft-clip baseQ sum + 13*hardclip_len
    int cap_mq;          // capped MQ on 0..C scale (internal)
    int cap_mq60;        // capped MQ rescaled to 0..60
    int eff_mq60;        // effective MQ = min(mapq_raw, cap_mq60)
    long long ins_len;   // total insertion length in read
    long long del_len;   // total deletion length in read
} ReadQC;

static void compute_read_qc(const bam1_t *b, int C, ReadQC *out){
    memset(out,0,sizeof(*out));
    out->mapq_raw = b->core.qual;

    const uint32_t *cig=bam_get_cigar(b); int nc=b->core.n_cigar;
    uint8_t *qual = bam_get_qual(b);
    int scL = left_clip_len(b,BAM_CSOFT_CLIP), scR = right_clip_len(b,BAM_CSOFT_CLIP);
    int hcL = left_clip_len(b,BAM_CHARD_CLIP), hcR = right_clip_len(b,BAM_CHARD_CLIP);
    int soft = scL + scR, hard = hcL + hcR;
    out->clipped_bases = soft + hard;

    // ClipQ
    long long sc_qsum=0;
    if (qual && scL){ for(int i=0;i<scL;i++) sc_qsum += qual[i]; }
    if (qual && scR){ int lq=b->core.l_qseq; for(int i=lq-scR;i<lq;i++) sc_qsum += qual[i]; }
    out->ClipQ = (int)(sc_qsum + 13LL*hard);

    // Count M_q13; tally indel lengths
    int rpos=0;
    for(int k=0;k<nc;k++){
        int op=bam_cigar_op(cig[k]), ln=bam_cigar_oplen(cig[k]);
        if (op==BAM_CMATCH||op==BAM_CEQUAL||op==BAM_CDIFF){
            for(int i=0;i<ln;i++){
                int q = qual ? qual[rpos+i] : 0;
                if (q >= BQ_MIN) out->M_q13++;
            }
            rpos += ln;
        } else if (op==BAM_CINS){ out->ins_len += ln; rpos += ln; }
        else if (op==BAM_CDEL){ out->del_len += ln; }
        else if (op==BAM_CSOFT_CLIP){ rpos += ln; }
    }

    // SubQ & X via MD (BQ>=13, cap 33)
    int mmn=0; Mismatch *mm = mismatches_from_MD(b, &mmn);
    for(int i=0;i<mmn;i++){
        if (mm[i].bq >= BQ_MIN){
            out->X_ge13++;
            out->SubQ += mm[i].bq < BQ_CAP ? mm[i].bq : BQ_CAP;
        }
    }
    if (mm) free(mm);

    // capped MQ (samtools -C math)
    // t = SubQ - 10*log10(M^X / X!) + ClipQ/5; Cap = sqrt((C - t)/C)*C  (or -1 if t > C)
    double t = 0.0;
    if (out->X_ge13 > 0 && out->M_q13 > 0){
        // log10(M^X / X!) = X*log10(M) - log10(X!)
        double log10_t = out->X_ge13 * log10((double)out->M_q13)
                       - ( log10(tgamma((double)out->X_ge13 + 1.0)) / log(10.0) );
        t = (double)out->SubQ - 10.0*log10_t + (double)out->ClipQ/5.0;
    } else {
        t = (double)out->ClipQ/5.0; // no subs
    }
    if (t > C) out->cap_mq = -1;
    else {
        if (t < 0) t = 0;
        double cap = sqrt((C - t)/C) * C;
        out->cap_mq = (int)(cap + 0.499);
    }
    // Rescale cap from 0..C to 0..60 (samtools-style)
    int cap_c = out->cap_mq < 0 ? 0 : out->cap_mq;     // treat -1 as 0 for rescaling
    int cap60 = (int)( (cap_c * 60.0) / (double)C + 0.5 );
    if (cap60 > 60) cap60 = 60;
    if (cap60 < 0)  cap60 = 0;
    out->cap_mq60 = cap60;

    // Effective MQ is usually limited by the cap
    out->eff_mq60 = (out->mapq_raw < out->cap_mq60) ? out->mapq_raw : out->cap_mq60;
}

// Apply one read to an interval [iv_start, iv_end):
//  - update per-site depth (count deletions as covering),
//  - push WR metrics once per site overlapped,
//  - add local INS/DEL length at start sites,
//  - return aligned_bp overlap (sum of M/= /X overlap lengths)
static long long apply_read_to_interval(const bam1_t *b, int iv_tid, int iv_start, int iv_end,
                                        const ReadQC *rq, Site *sites)
{
    const uint32_t *cig=bam_get_cigar(b); int nc=b->core.n_cigar;
    if (b->core.tid != iv_tid) return 0;
    int64_t ref = b->core.pos; // 0-based
    int rpos = 0;
    long long aligned_bp = 0;

    for(int k=0;k<nc;k++){
        int op=bam_cigar_op(cig[k]), ln=bam_cigar_oplen(cig[k]);
        if (op==BAM_CMATCH||op==BAM_CEQUAL||op==BAM_CDIFF){
            // per-base over ref
            for(int i=0;i<ln;i++){
                int64_t rp = ref + i;
                if (rp >= iv_start && rp < iv_end){
                    int idx = (int)(rp - iv_start);
                    // depth++
                    sites[idx].depth += 1;
                    // WR metrics once per site
                    vpush(&sites[idx].mq,    (double)rq->mapq_raw);
                    vpush(&sites[idx].capmq, (double)rq->cap_mq60);
                    vpush(&sites[idx].effmq, (double)rq->eff_mq60);
                    vpush(&sites[idx].subq,  (double)rq->SubQ);
                    vpush(&sites[idx].clipq, (double)rq->ClipQ);
                    sites[idx].clip_bases_sum += rq->clipped_bases;
                }
            }
            // aligned overlap length for the interval
            int64_t segL = ref, segR = ref + ln; // [segL, segR)
            int64_t ovL = segL > iv_start ? segL : iv_start;
            int64_t ovR = segR < iv_end   ? segR : iv_end;
            if (ovR > ovL) aligned_bp += (ovR - ovL);
            ref  += ln; rpos += ln;
        } else if (op==BAM_CDEL){
            // deletions cover reference positions for depth
            for(int i=0;i<ln;i++){
                int64_t rp = ref + i;
                if (rp >= iv_start && rp < iv_end){
                    int idx = (int)(rp - iv_start);
                    sites[idx].depth += 1;
                }
            }
            // assign local DEL length at start site
            if (ref >= iv_start && ref < iv_end){
                int idx = (int)(ref - iv_start);
                sites[idx].del_len_sum += ln;
            }
            ref += ln;
        } else if (op==BAM_CINS){
            // assign local INS length at current site (before ref advance)
            if (ref >= iv_start && ref < iv_end){
                int idx = (int)(ref - iv_start);
                sites[idx].ins_len_sum += ln;
            }
            rpos += ln; // insertion consumes read only
        } else if (op==BAM_CREF_SKIP){ ref += ln; }
        else if (op==BAM_CSOFT_CLIP){ rpos += ln; }
        else if (op==BAM_CHARD_CLIP){ /* nothing */ }
        else if (op==BAM_CPAD){ /* nothing */ }
    }
    return aligned_bp;
}

// Add local mismatches (from MD) to per-site counters
static void add_local_mismatches(const bam1_t *b, int iv_tid, int iv_start, int iv_end, Site *sites){
    if (b->core.tid != iv_tid) return;
    // build array of ref positions corresponding to aligned read positions (M/= /X)
    const uint32_t *cig=bam_get_cigar(b); int nc=b->core.n_cigar;
    int64_t ref = b->core.pos; int rpos=0, cap=0;
    for(int k=0;k<nc;k++){
        int op=bam_cigar_op(cig[k]), ln=bam_cigar_oplen(cig[k]);
        if (op==BAM_CMATCH||op==BAM_CEQUAL||op==BAM_CDIFF) cap += ln;
        else if (op==BAM_CINS||op==BAM_CSOFT_CLIP) rpos += ln;
        else if (op==BAM_CDEL||op==BAM_CREF_SKIP) ref += ln;
    }
    int64_t *refpos = cap? (int64_t*)malloc(cap*sizeof(int64_t)) : NULL;
    if (cap && !refpos){ perror("malloc"); exit(1); }
    ref = b->core.pos; rpos=0; int idx=0;
    for(int k=0;k<nc;k++){
        int op=bam_cigar_op(cig[k]), ln=bam_cigar_oplen(cig[k]);
        if (op==BAM_CMATCH||op==BAM_CEQUAL||op==BAM_CDIFF){
            for(int i=0;i<ln;i++) refpos[idx++] = ref + i;
            ref += ln; rpos += ln;
        } else if (op==BAM_CINS){ rpos += ln; }
        else if (op==BAM_CDEL||op==BAM_CREF_SKIP){ ref += ln; }
        else if (op==BAM_CSOFT_CLIP){ rpos += ln; }
    }
    // Parse MD to mark mismatches
    uint8_t *mdp = bam_aux_get(b,"MD");
    if (mdp){
        const char *md = bam_aux2Z(mdp);
        int ai=0; const char *p=md;
        while (*p){
            if (isdigit((unsigned char)*p)){
                int v=0; while (isdigit((unsigned char)*p)){ v = v*10 + (*p-'0'); ++p; }
                ai += v;
            } else if (*p=='^'){ ++p; while (*p && isalpha((unsigned char)*p)) ++p; }
            else if (isalpha((unsigned char)*p)){
                if (ai < cap){
                    int64_t rp = refpos[ai];
                    if (rp >= iv_start && rp < iv_end){
                        int sidx = (int)(rp - iv_start);
                        sites[sidx].mismatches += 1;
                    }
                    ai++;
                }
                ++p;
            } else ++p;
        }
    }
    if (refpos) free(refpos);
}

// ---------- BED reader (BED3+ minimal) ----------
static BedIv *read_bed(const char *path, int *n_out){
    FILE *f = strcmp(path,"-")==0? stdin : fopen(path,"r");
    if(!f){ fprintf(stderr,"ERROR: cannot open BED %s\n", path); exit(1); }
    int cap=1024, n=0;
    BedIv *a = (BedIv*)malloc(cap*sizeof(BedIv));
    if(!a){ perror("malloc"); exit(1); }
    char buf[1<<16];
    while (fgets(buf,sizeof(buf),f)){
        if (buf[0]=='#' || buf[0]=='\n') continue;
        char chr[1024]; long s,e;
        if (sscanf(buf, "%1023s %ld %ld", chr, &s, &e) != 3) continue;
        if (n==cap){ cap*=2; a=(BedIv*)realloc(a,cap*sizeof(BedIv)); if(!a){ perror("realloc"); exit(1); } }
        a[n].chr = xstrdup(chr);
        a[n].start = (int)s; a[n].end = (int)e; n++;
    }
    if (f!=stdin) fclose(f);
    *n_out = n; return a;
}
static void free_bed(BedIv *a, int n){ for(int i=0;i<n;i++) free(a[i].chr); free(a); }

// ---------- MAIN ----------
static void usage(void){
    fprintf(stderr,
      "Usage: mq_forensics -b in.bam -r regions.bed -C INT -d MIN_DEPTH -o per_site.tsv -O per_interval.tsv\n");
}

int main(int argc, char **argv){
    const char *bam_path=NULL, *bed_path=NULL, *out_site=NULL, *out_iv=NULL;
    int C = 50;
    int min_depth = 0; // per-site minimum depth threshold

    for (int i=1;i<argc;i++){
        if (!strcmp(argv[i],"-b") && i+1<argc) bam_path=argv[++i];
        else if (!strcmp(argv[i],"-r") && i+1<argc) bed_path=argv[++i];
        else if (!strcmp(argv[i],"-C") && i+1<argc) C = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-d") && i+1<argc) min_depth = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-o") && i+1<argc) out_site=argv[++i];
        else if (!strcmp(argv[i],"-O") && i+1<argc) out_iv=argv[++i];
        else { usage(); return 2; }
    }
    if (!bam_path || !bed_path || !out_site || !out_iv){ usage(); return 2; }

    samFile *fp = sam_open(bam_path,"r");
    if(!fp){ fprintf(stderr,"ERROR: open BAM %s\n", bam_path); return 1; }
    bam_hdr_t *hdr = sam_hdr_read(fp);
    if(!hdr){ fprintf(stderr,"ERROR: read header\n"); sam_close(fp); return 1; }
    hts_idx_t *idx = sam_index_load(fp, bam_path);
    if(!idx){ fprintf(stderr,"ERROR: BAM index not found for %s\n", bam_path); return 1; }

    int nbed=0; BedIv *bed = read_bed(bed_path, &nbed);
    FILE *fo_site = fopen(out_site,"w");
    FILE *fo_iv   = fopen(out_iv,"w");
    if(!fo_site||!fo_iv){ fprintf(stderr,"ERROR: open outputs\n"); return 1; }

    // headers
    fprintf(fo_site, "chrom\tpos\tdepth\tmismatch_bases\tins_len_sum\tdel_len_sum\t"
                     "clip_bases_sum\t"
                     "mq_mean\tmq_median\tmq_sd\t"
                     "capmq60_mean\tcapmq60_median\tcapmq60_sd\t"
                     "effmq_mean\teffmq_median\teffmq_sd\t"
                     "subQ_mean\tsubQ_median\tsubQ_sd\t"
                     "clipQ_mean\tclipQ_median\tclipQ_sd\n");
    fprintf(fo_iv,   "chrom\tstart\tend\taligned_bp_in_interval\n");

    bam1_t *b = bam_init1();

    for (int iv=0; iv<nbed; iv++){
        int tid = sam_hdr_name2tid(hdr, bed[iv].chr);
        if (tid < 0){
            fprintf(stderr,"WARN: contig %s not in BAM header; skipping %s:%d-%d\n",
                    bed[iv].chr, bed[iv].chr, bed[iv].start, bed[iv].end);
            continue;
        }
        int L = bed[iv].end - bed[iv].start;
        if (L <= 0) continue;

        // alloc per-site
        Site *sites = (Site*)calloc(L, sizeof(Site));
        if(!sites){ perror("calloc"); exit(1); }

        // iterate reads overlapping interval
        hts_itr_t *itr = sam_itr_queryi(idx, tid, bed[iv].start, bed[iv].end);
        if (!itr){
            fprintf(stderr,"WARN: bad iterator for %s:%d-%d\n", bed[iv].chr, bed[iv].start, bed[iv].end);
            free(sites); continue;
        }
        long long aligned_bp_sum = 0;

        while (sam_itr_next(fp, itr, b) >= 0){
            if (b->core.flag & (BAM_FUNMAP|BAM_FSECONDARY|BAM_FSUPPLEMENTARY)) continue;

            // compute per-read metrics once
            ReadQC rq; compute_read_qc(b, C, &rq);

            // apply to interval: depth, WR per-site, INS/DEL sums, aligned_bp
            aligned_bp_sum += apply_read_to_interval(b, tid, bed[iv].start, bed[iv].end, &rq, sites);

            // local mismatches per site (from MD)
            add_local_mismatches(b, tid, bed[iv].start, bed[iv].end, sites);
        }
        hts_itr_destroy(itr);

        // write per-site rows
        for (int i=0;i<L;i++){
            int pos1 = bed[iv].start + i + 1; // 1-based

            if (sites[i].depth < min_depth) {
                // depth below threshold → print NA for everything numeric
                // (there are 20 numeric columns after pos)
                fprintf(fo_site,
                    "%s\t%d\t"
                    "NA\tNA\tNA\tNA\tNA\t"
                    "NA\tNA\tNA\t"
                    "NA\tNA\tNA\t"
                    "NA\tNA\tNA\t"
                    "NA\tNA\tNA\t"
                    "NA\tNA\tNA\n",
                    bed[iv].chr, pos1);
            } else {
                double mq_mu   = dmean(&sites[i].mq);
                double mq_med  = dmedian(&sites[i].mq);
                double mq_sdv  = dsd(&sites[i].mq);

                double cmu     = dmean(&sites[i].capmq);
                double cmed    = dmedian(&sites[i].capmq);
                double csd     = dsd(&sites[i].capmq);

                double emu     = dmean(&sites[i].effmq);
                double emed    = dmedian(&sites[i].effmq);
                double esd     = dsd(&sites[i].effmq);

                double sq_mu   = dmean(&sites[i].subq);
                double sq_med  = dmedian(&sites[i].subq);
                double sq_sdv  = dsd(&sites[i].subq);

                double cq_mu   = dmean(&sites[i].clipq);
                double cq_med  = dmedian(&sites[i].clipq);
                double cq_sdv  = dsd(&sites[i].clipq);

                fprintf(fo_site,
                    "%s\t%d\t%d\t%d\t%lld\t%lld\t%lld\t"
                    "%.3f\t%.3f\t%.3f\t"
                    "%.3f\t%.3f\t%.3f\t"
                    "%.3f\t%.3f\t%.3f\t"
                    "%.3f\t%.3f\t%.3f\t"
                    "%.3f\t%.3f\t%.3f\n",
                    bed[iv].chr, pos1,
                    sites[i].depth, sites[i].mismatches,
                    sites[i].ins_len_sum, sites[i].del_len_sum, sites[i].clip_bases_sum,
                    mq_mu, mq_med, mq_sdv,
                    cmu,   cmed,   csd,
                    emu,   emed,   esd,
                    sq_mu, sq_med, sq_sdv,
                    cq_mu, cq_med, cq_sdv);
            }
            // free per-site vectors
            free(sites[i].mq.v);
            free(sites[i].capmq.v);
            free(sites[i].effmq.v);
            free(sites[i].subq.v);
            free(sites[i].clipq.v);
        }
        // per-interval row
        fprintf(fo_iv, "%s\t%d\t%d\t%lld\n", bed[iv].chr, bed[iv].start, bed[iv].end, aligned_bp_sum);

        free(sites);
    }

    bam_destroy1(b);
    fclose(fo_site); fclose(fo_iv);
    free_bed(bed, nbed);
    hts_idx_destroy(idx);
    bam_hdr_destroy(hdr);
    sam_close(fp);
    return 0;
}
