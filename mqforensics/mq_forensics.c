// mq_forensics.c — combined tool
// Modes:
//   1) stats (default): BAM+BED per-site & per-interval stats (WR attribution, -C cap, capmq60, effmq, NA with -d)
//   2) summarize: streaming site-level summary across samples (sorted TSV), producing depth-normalized rates + direct metrics.
//
// Build:
//   gcc -O3 -std=c11 mq_forensics.c -o mq_forensics $(pkg-config --cflags --libs htslib) -lm
//
// Usage (stats):
//   mq_forensics -b in.bam -r regions.bed -C INT -d MIN_DEPTH -o per_site.tsv -O per_interval.tsv
//
// Usage (summarize; input must be sorted by chrom,pos):
//   mq_forensics summarize -i all_samples.sorted.tsv -o per_site_across_samples.tsv
//   # or stream: ... | mq_forensics summarize -i - -o -
//
// Notes:
//   - Requires MD/NM in BAM for mismatch parsing. If absent: samtools calmd -bAr ref.fa in.bam > out.bam
//   - WR attribution: each overlapping read contributes its {MAPQ, capmq60, effmq, SubQ, ClipQ, clipped_bases} to every site it overlaps.
//   - capmq60 is the -C cap rescaled to 0..60; effmq = min(raw MAPQ, capmq60).

#define _GNU_SOURCE
#include <htslib/sam.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/*========================== Common helpers ==========================*/

#define BQ_MIN 13
#define BQ_CAP 33

static char* xstrdup(const char* s){
    size_t n = strlen(s) + 1;
    char *p = (char*)malloc(n);
    if (!p) { perror("malloc"); exit(1); }
    memcpy(p, s, n);
    return p;
}

typedef struct { double *v; int n, cap; } vecd;
static inline void vpush(vecd *x, double a){
    if(x->n==x->cap){
        x->cap = x->cap? x->cap*2:256;
        x->v = (double*)realloc(x->v, x->cap*sizeof(double));
        if(!x->v){ perror("realloc"); exit(1); }
    }
    x->v[x->n++]=a;
}
static int dcmp(const void*a,const void*b){
    double x=*(const double*)a, y=*(const double*)b;
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

/*========================== stats mode (BAM→per-site) ==========================*/

typedef struct {
    int rp; int bq;
} Mismatch;

static Mismatch* mismatches_from_MD(const bam1_t *b, int *n_out){
    uint8_t *mdp = bam_aux_get(b,"MD");
    if (!mdp){ *n_out=0; return NULL; }
    const char *md = bam_aux2Z(mdp); if(!md){ *n_out=0; return NULL; }

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
    Mismatch *mm = cap? (Mismatch*)malloc(cap*sizeof(Mismatch)) : NULL;
    if (cap && !mm){ perror("malloc"); exit(1); }
    int mmn=0;

    while (*p){
        if (isdigit((unsigned char)*p)){
            int v=0; while (isdigit((unsigned char)*p)){ v = v*10 + (*p-'0'); ++p; }
            ai += v;
        } else if (*p=='^'){
            ++p; while (*p && isalpha((unsigned char)*p)) ++p;
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

static inline int left_clip_len(const bam1_t *b, int op){
    const uint32_t *cig=bam_get_cigar(b);
    if (b->core.n_cigar==0) return 0;
    return (bam_cigar_op(cig[0])==op) ? bam_cigar_oplen(cig[0]) : 0;
}
static inline int right_clip_len(const bam1_t *b, int op){
    const uint32_t *cig=bam_get_cigar(b); int n=b->core.n_cigar; if(!n) return 0;
    return (bam_cigar_op(cig[n-1])==op) ? bam_cigar_oplen(cig[n-1]) : 0;
}

typedef struct {
    int mapq_raw;
    int clipped_bases;
    int M_q13;
    int X_ge13;
    int SubQ;
    int ClipQ;
    int cap_mq;     // 0..C or -1 if t > C
    int cap_mq60;   // rescaled 0..60 (for reporting)
    int eff_mq60;   // min(raw, cap_mq60)
    long long ins_len;
    long long del_len;
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

    long long sc_qsum=0;
    if (qual && scL){ for(int i=0;i<scL;i++) sc_qsum += qual[i]; }
    if (qual && scR){ int lq=b->core.l_qseq; for(int i=lq-scR;i<lq;i++) sc_qsum += qual[i]; }
    out->ClipQ = (int)(sc_qsum + 13LL*hard);

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

    int mmn=0; Mismatch *mm = mismatches_from_MD(b, &mmn);
    for(int i=0;i<mmn;i++){
        if (mm[i].bq >= BQ_MIN){
            out->X_ge13++;
            out->SubQ += mm[i].bq < BQ_CAP ? mm[i].bq : BQ_CAP;
        }
    }
    if (mm) free(mm);

    double t = 0.0;
    if (out->X_ge13 > 0 && out->M_q13 > 0){
        double log10_t = out->X_ge13 * log10((double)out->M_q13)
                       - ( log10(tgamma((double)out->X_ge13 + 1.0)) / log(10.0) );
        t = (double)out->SubQ - 10.0*log10_t + (double)out->ClipQ/5.0;
    } else {
        t = (double)out->ClipQ/5.0;
    }
    if (t > C) out->cap_mq = -1;
    else {
        if (t < 0) t = 0;
        double cap = sqrt((C - t)/C) * C;
        out->cap_mq = (int)(cap + 0.499);
    }
    int cap_c = out->cap_mq < 0 ? 0 : out->cap_mq;
    int cap60 = (int)( (cap_c * 60.0) / (double)C + 0.5 );
    if (cap60 > 60) cap60 = 60;
    if (cap60 < 0)  cap60 = 0;
    out->cap_mq60 = cap60;
    out->eff_mq60 = (out->mapq_raw < out->cap_mq60) ? out->mapq_raw : out->cap_mq60;
}

typedef struct {
    int depth;
    int mismatches;
    long long ins_len_sum;
    long long del_len_sum;
    long long clip_bases_sum;
    vecd mq;
    vecd capmq;
    vecd effmq;
    vecd subq;
    vecd clipq;
} Site;

typedef struct { char *chr; int start, end; } BedIv;

static long long apply_read_to_interval(const bam1_t *b, int iv_tid, int iv_start, int iv_end,
                                        const ReadQC *rq, Site *sites)
{
    const uint32_t *cig=bam_get_cigar(b); int nc=b->core.n_cigar;
    if (b->core.tid != iv_tid) return 0;
    int64_t ref = b->core.pos;
    int rpos = 0;
    long long aligned_bp = 0;

    for(int k=0;k<nc;k++){
        int op=bam_cigar_op(cig[k]), ln=bam_cigar_oplen(cig[k]);
        if (op==BAM_CMATCH||op==BAM_CEQUAL||op==BAM_CDIFF){
            for(int i=0;i<ln;i++){
                int64_t rp = ref + i;
                if (rp >= iv_start && rp < iv_end){
                    int idx = (int)(rp - iv_start);
                    sites[idx].depth += 1;
                    vpush(&sites[idx].mq,    (double)rq->mapq_raw);
                    vpush(&sites[idx].capmq, (double)rq->cap_mq60);
                    vpush(&sites[idx].effmq, (double)rq->eff_mq60);
                    vpush(&sites[idx].subq,  (double)rq->SubQ);
                    vpush(&sites[idx].clipq, (double)rq->ClipQ);
                    sites[idx].clip_bases_sum += rq->clipped_bases;
                }
            }
            int64_t segL = ref, segR = ref + ln;
            int64_t ovL = segL > iv_start ? segL : iv_start;
            int64_t ovR = segR < iv_end   ? segR : iv_end;
            if (ovR > ovL) aligned_bp += (ovR - ovL);
            ref  += ln; rpos += ln;
        } else if (op==BAM_CDEL){
            for(int i=0;i<ln;i++){
                int64_t rp = ref + i;
                if (rp >= iv_start && rp < iv_end){
                    int idx = (int)(rp - iv_start);
                    sites[idx].depth += 1;
                }
            }
            if (ref >= iv_start && ref < iv_end){
                int idx = (int)(ref - iv_start);
                sites[idx].del_len_sum += ln;
            }
            ref += ln;
        } else if (op==BAM_CINS){
            if (ref >= iv_start && ref < iv_end){
                int idx = (int)(ref - iv_start);
                sites[idx].ins_len_sum += ln;
            }
            rpos += ln;
        } else if (op==BAM_CREF_SKIP){ ref += ln; }
        else if (op==BAM_CSOFT_CLIP){ rpos += ln; }
    }
    return aligned_bp;
}

static void add_local_mismatches(const bam1_t *b, int iv_tid, int iv_start, int iv_end, Site *sites){
    if (b->core.tid != iv_tid) return;
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

/*========================== BED reader ==========================*/

typedef struct { BedIv *a; int n, cap; } BedVec;
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

/*========================== summarize mode (TSV→per-site across samples) ==========================*/

static int split_tabs(char *line, char **out, int max){
    int n=0; char *p=line, *start=p;
    while (*p){
        if (*p=='\t' || *p=='\n'){
            *p=0; out[n++]=start; start=p+1;
            if (n==max) break;
        }
        p++;
    }
    if (n<max && *start) out[n++]=start;
    return n;
}
static int isNA(const char *s){ return s[0]=='N' && s[1]=='A' && s[2]==0; }

/* Direct metrics summarized as-is if present */
static const char *DIRECT_COLS[] = {
  "mq_mean","mq_median","mq_sd",
  "capmq60_mean","capmq60_median","capmq60_sd",
  "effmq_mean","effmq_median","effmq_sd",
  "subQ_mean","subQ_median","subQ_sd",
  "clipQ_mean","clipQ_median","clipQ_sd"
};
static const int NDIRECT = (int)(sizeof(DIRECT_COLS)/sizeof(DIRECT_COLS[0]));

typedef struct { const char *out_name; const char *sum_col; int sum_idx; vecd vals; } RateMetric;
static RateMetric RATE_COLS[] = {
  { "mismatch_rate", "mismatch_bases",  -1, {NULL,0,0} },
  { "ins_rate",      "ins_len_sum",     -1, {NULL,0,0} },
  { "del_rate",      "del_len_sum",     -1, {NULL,0,0} },
  { "clip_rate",     "clip_bases_sum",  -1, {NULL,0,0} }
};
static const int NRATE = (int)(sizeof(RATE_COLS)/sizeof(RATE_COLS[0]));

typedef struct {
    char name[64];
    int index;
    vecd vals;
} DirectMetric;

static int summarize_sites_stream(const char *in_path, const char *out_path){
    FILE *fi = stdin;
    if (in_path && strcmp(in_path,"-")!=0){
        fi = fopen(in_path,"r");
        if (!fi){ fprintf(stderr,"ERROR: open %s\n", in_path); return 1; }
    }
    FILE *fo = stdout;
    if (out_path && strcmp(out_path,"-")!=0){
        fo = fopen(out_path,"w");
        if (!fo){ fprintf(stderr,"ERROR: open %s\n", out_path); if (fi!=stdin) fclose(fi); return 1; }
    }

    char *line=NULL; size_t cap=0; ssize_t nread = getline(&line,&cap,fi);
    if (nread<0){ fprintf(stderr,"ERROR: empty input\n"); goto fail; }

    char *cols[1024]; int ncols = split_tabs(line, cols, 1024);
    if (ncols < 3){ fprintf(stderr,"ERROR: need at least chrom,pos and metrics\n"); goto fail; }

    int chrom_idx=0, pos_idx=1, depth_idx=-1;
    for (int j=0;j<ncols;j++){ if (!strcmp(cols[j],"depth")) { depth_idx=j; break; } }
    if (depth_idx<0){ fprintf(stderr,"ERROR: 'depth' column not found (needed for rates)\n"); goto fail; }

    DirectMetric DM[128]; int K=0;
    for (int i=0;i<NDIRECT;i++){
        for (int j=0;j<ncols;j++){
            if (!strcmp(DIRECT_COLS[i], cols[j])){
                strncpy(DM[K].name, DIRECT_COLS[i], sizeof(DM[K].name)-1);
                DM[K].name[sizeof(DM[K].name)-1]=0; DM[K].index=j;
                DM[K].vals.v=NULL; DM[K].vals.n=DM[K].vals.cap=0; K++; break;
            }
        }
    }
    int any_rate=0;
    for (int r=0;r<NRATE;r++){
        RATE_COLS[r].sum_idx = -1;
        for (int j=0;j<ncols;j++){
            if (!strcmp(RATE_COLS[r].sum_col, cols[j])){
                RATE_COLS[r].sum_idx=j;
                RATE_COLS[r].vals.v=NULL; RATE_COLS[r].vals.n=RATE_COLS[r].vals.cap=0;
                any_rate=1; break;
            }
        }
    }
    if (!any_rate) fprintf(stderr,"WARN: no sum columns for rate metrics present; only direct metrics will be summarized\n");

    // header
    fprintf(fo,"chrom\tpos");
    for (int r=0;r<NRATE;r++){
        if (RATE_COLS[r].sum_idx>=0)
            fprintf(fo,"\t%s_mean\t%s_median\t%s_sd", RATE_COLS[r].out_name, RATE_COLS[r].out_name, RATE_COLS[r].out_name);
    }
    for (int k=0;k<K;k++){
        fprintf(fo,"\t%s_mean\t%s_median\t%s_sd", DM[k].name, DM[k].name, DM[k].name);
    }
    fprintf(fo,"\n");

    char *fields[2048];
    char prev_chr[512]=""; long prev_pos=-1;

    while ((nread=getline(&line,&cap,fi))>=0){
        if (line[0]=='\n') continue;
        int nf = split_tabs(line, fields, 2048);
        if (nf < ncols) continue;

        char *chr = fields[chrom_idx];
        long pos = strtol(fields[pos_idx], NULL, 10);

        if (prev_pos!=-1 && (pos!=prev_pos || strcmp(chr,prev_chr)!=0)){
            fprintf(fo,"%s\t%ld", prev_chr, prev_pos);
            for (int r=0;r<NRATE;r++){
                if (RATE_COLS[r].sum_idx>=0){
                    if (RATE_COLS[r].vals.n==0) fprintf(fo,"\tNA\tNA\tNA");
                    else {
                        double mu=dmean(&RATE_COLS[r].vals), md=dmedian(&RATE_COLS[r].vals), sdv=dsd(&RATE_COLS[r].vals);
                        fprintf(fo,"\t%.6f\t%.6f\t%.6f", mu, md, sdv);
                    }
                    free(RATE_COLS[r].vals.v); RATE_COLS[r].vals.v=NULL; RATE_COLS[r].vals.n=RATE_COLS[r].vals.cap=0;
                }
            }
            for (int k=0;k<K;k++){
                if (DM[k].vals.n==0) fprintf(fo,"\tNA\tNA\tNA");
                else {
                    double mu=dmean(&DM[k].vals), md=dmedian(&DM[k].vals), sdv=dsd(&DM[k].vals);
                    fprintf(fo,"\t%.6f\t%.6f\t%.6f", mu, md, sdv);
                }
                free(DM[k].vals.v); DM[k].vals.v=NULL; DM[k].vals.n=DM[k].vals.cap=0;
            }
            fprintf(fo,"\n");
        }
        if (prev_pos==-1 || pos!=prev_pos || strcmp(chr,prev_chr)!=0){
            strncpy(prev_chr, chr, sizeof(prev_chr)-1); prev_chr[sizeof(prev_chr)-1]=0; prev_pos=pos;
        }

        // rates: sum/depth per row
        double depth=0.0; int depth_ok=0;
        if (!isNA(fields[depth_idx])){
            char *e=NULL; depth=strtod(fields[depth_idx], &e);
            if (e!=fields[depth_idx] && depth>0.0) depth_ok=1;
        }
        if (depth_ok){
            for (int r=0;r<NRATE;r++){
                int si=RATE_COLS[r].sum_idx;
                if (si>=0 && !isNA(fields[si])){
                    char *e=NULL; double sumv=strtod(fields[si], &e);
                    if (e!=fields[si]) vpush(&RATE_COLS[r].vals, sumv/depth);
                }
            }
        }
        // direct metrics
        for (int k=0;k<K;k++){
            const char *s = fields[DM[k].index];
            if (!isNA(s)){
                char *e=NULL; double v=strtod(s,&e);
                if (e!=s) vpush(&DM[k].vals, v);
            }
        }
    }
    if (prev_pos!=-1){
        fprintf(fo,"%s\t%ld", prev_chr, prev_pos);
        for (int r=0;r<NRATE;r++){
            if (RATE_COLS[r].sum_idx>=0){
                if (RATE_COLS[r].vals.n==0) fprintf(fo,"\tNA\tNA\tNA");
                else {
                    double mu=dmean(&RATE_COLS[r].vals), md=dmedian(&RATE_COLS[r].vals), sdv=dsd(&RATE_COLS[r].vals);
                    fprintf(fo,"\t%.6f\t%.6f\t%.6f", mu, md, sdv);
                }
                free(RATE_COLS[r].vals.v);
            }
        }
        for (int k=0;k<K;k++){
            if (DM[k].vals.n==0) fprintf(fo,"\tNA\tNA\tNA");
            else {
                double mu=dmean(&DM[k].vals), md=dmedian(&DM[k].vals), sdv=dsd(&DM[k].vals);
                fprintf(fo,"\t%.6f\t%.6f\t%.6f", mu, md, sdv);
            }
            free(DM[k].vals.v);
        }
        fprintf(fo,"\n");
    }

    free(line);
    if (fo!=stdout) fclose(fo);
    if (fi!=stdin) fclose(fi);
    return 0;

fail:
    free(line);
    return 2;
}

/*========================== usage & main ==========================*/

static void usage(void){
    fprintf(stderr,
      "Usage (stats mode):\n"
      "  mq_forensics -b in.bam -r regions.bed -C INT -d MIN_DEPTH -o per_site.tsv -O per_interval.tsv\n"
      "\n"
      "Usage (summarize mode; input sorted by chrom,pos):\n"
      "  mq_forensics summarize -i all_samples.sorted.tsv -o per_site_across_samples.tsv\n"
      "  (use '-' for stdin/stdout)\n");
}

int main(int argc, char **argv){
    // subcommand: summarize
    if (argc>1 && strcmp(argv[1],"summarize")==0){
        const char *in_path="-", *out_path="-";
        for (int i=2;i<argc;i++){
            if (!strcmp(argv[i],"-i") && i+1<argc) in_path=argv[++i];
            else if (!strcmp(argv[i],"-o") && i+1<argc) out_path=argv[++i];
            else { usage(); return 2; }
        }
        return summarize_sites_stream(in_path, out_path);
    }

    const char *bam_path=NULL, *bed_path=NULL, *out_site=NULL, *out_iv=NULL;
    int C = 50;
    int min_depth = 0;

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

        Site *sites = (Site*)calloc(L, sizeof(Site));
        if(!sites){ perror("calloc"); exit(1); }

        hts_itr_t *itr = sam_itr_queryi(idx, tid, bed[iv].start, bed[iv].end);
        if (!itr){
            fprintf(stderr,"WARN: bad iterator for %s:%d-%d\n", bed[iv].chr, bed[iv].start, bed[iv].end);
            free(sites); continue;
        }
        long long aligned_bp_sum = 0;

        while (sam_itr_next(fp, itr, b) >= 0){
            if (b->core.flag & (BAM_FUNMAP|BAM_FSECONDARY|BAM_FSUPPLEMENTARY)) continue;
            ReadQC rq; compute_read_qc(b, C, &rq);
            aligned_bp_sum += apply_read_to_interval(b, tid, bed[iv].start, bed[iv].end, &rq, sites);
            add_local_mismatches(b, tid, bed[iv].start, bed[iv].end, sites);
        }
        hts_itr_destroy(itr);

        for (int i=0;i<L;i++){
            int pos1 = bed[iv].start + i + 1;
            if (sites[i].depth < min_depth){
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
                double mq_mu=dmean(&sites[i].mq), mq_med=dmedian(&sites[i].mq), mq_sdv=dsd(&sites[i].mq);
                double cmu=dmean(&sites[i].capmq), cmed=dmedian(&sites[i].capmq), csd=dsd(&sites[i].capmq);
                double emu=dmean(&sites[i].effmq), emed=dmedian(&sites[i].effmq), esd=dsd(&sites[i].effmq);
                double sq_mu=dmean(&sites[i].subq), sq_med=dmedian(&sites[i].subq), sq_sdv=dsd(&sites[i].subq);
                double cq_mu=dmean(&sites[i].clipq), cq_med=dmedian(&sites[i].clipq), cq_sdv=dsd(&sites[i].clipq);
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
            free(sites[i].mq.v);
            free(sites[i].capmq.v);
            free(sites[i].effmq.v);
            free(sites[i].subq.v);
            free(sites[i].clipq.v);
        }
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
