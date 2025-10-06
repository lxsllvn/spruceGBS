#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <htslib/sam.h>
#include "common.h"

// xstrdup
static char* xstrdup(const char* s){
    size_t n=strlen(s)+1; char *p=(char*)malloc(n);
    if(!p){ perror("malloc"); exit(1); }
    memcpy(p,s,n); return p;
}

// BED reader
int read_bed(const char *path, BedIv **out){
    FILE *f = strcmp(path,"-")==0? stdin : fopen(path,"r");
    if(!f){ fprintf(stderr,"ERROR: cannot open BED %s\n", path); exit(1); }
    int cap=1024,n=0; BedIv *a=(BedIv*)malloc(cap*sizeof(BedIv));
    if(!a){ perror("malloc"); exit(1); }
    char buf[1<<16];
    while (fgets(buf,sizeof(buf),f)){
        if (buf[0]=='#' || buf[0]=='\n') continue;
        char chr[1024]; long s,e;
        if (sscanf(buf,"%1023s %ld %ld", chr, &s, &e) != 3) continue;
        if (n==cap){ cap*=2; a=(BedIv*)realloc(a,cap*sizeof(BedIv)); if(!a){ perror("realloc"); exit(1);} }
        a[n].chr=xstrdup(chr); a[n].start=(int)s; a[n].end=(int)e; n++;
    }
    if (f!=stdin) fclose(f);
    *out=a; return n;
}
void free_bed(BedIv *a, int n){ for(int i=0;i<n;i++) free(a[i].chr); free(a); }

// helpers
static inline int left_clip_len(const bam1_t *b, int op){
    const uint32_t *cig=bam_get_cigar(b);
    if (b->core.n_cigar==0) return 0;
    return ((int)bam_cigar_op(cig[0]) == op) ? (int)bam_cigar_oplen(cig[0]) : 0;
}
static inline int right_clip_len(const bam1_t *b, int op){
    const uint32_t *cig=bam_get_cigar(b); int n=b->core.n_cigar; if(!n) return 0;
    return ((int)bam_cigar_op(cig[n-1]) == op) ? (int)bam_cigar_oplen(cig[n-1]) : 0;
}

// mismatches via MD
typedef struct { int rp; int bq; } Mismatch;

static Mismatch* mismatches_from_MD(const bam1_t *b, int *n_out){
    uint8_t *mdp=bam_aux_get(b,"MD"); if(!mdp){ *n_out=0; return NULL; }
    const char *md=bam_aux2Z(mdp); if(!md){ *n_out=0; return NULL; }

    const uint32_t *cig=bam_get_cigar(b);
    int nc=b->core.n_cigar, rpos=0, cap=0;
    for(int k=0;k<nc;k++){
        int op=bam_cigar_op(cig[k]), ln=bam_cigar_oplen(cig[k]);
        if (op==BAM_CMATCH||op==BAM_CEQUAL||op==BAM_CDIFF) cap+=ln;
        else if (op==BAM_CINS||op==BAM_CSOFT_CLIP) rpos+=ln;
    }
    int *aligned = cap? (int*)malloc(cap*sizeof(int)) : NULL;
    if (cap && !aligned){ perror("malloc"); exit(1); }
    int idx=0; rpos=0;
    for(int k=0;k<nc;k++){
        int op=bam_cigar_op(cig[k]), ln=bam_cigar_oplen(cig[k]);
        if (op==BAM_CMATCH||op==BAM_CEQUAL||op==BAM_CDIFF){ for(int i=0;i<ln;i++) aligned[idx++]=rpos++; }
        else if (op==BAM_CINS||op==BAM_CSOFT_CLIP){ rpos+=ln; }
    }
    uint8_t *qual=bam_get_qual(b);
    Mismatch *mm = cap? (Mismatch*)malloc(cap*sizeof(Mismatch)) : NULL; if(cap && !mm){ perror("malloc"); exit(1); }
    int mmn=0; const char *p=md; int ai=0;
    while (*p){
        if (isdigit((unsigned char)*p)){
            int v=0; while(isdigit((unsigned char)*p)){ v=v*10+(*p-'0'); ++p; } ai+=v;
        } else if (*p=='^'){ ++p; while(*p && isalpha((unsigned char)*p)) ++p; }
        else if (isalpha((unsigned char)*p)){
            if (ai<cap){
                int rp=aligned[ai], bq=qual?qual[rp]:0;
                mm[mmn].rp=rp; mm[mmn].bq=bq; mmn++; ai++;
            }
            ++p;
        } else ++p;
    }
    free(aligned); *n_out=mmn; return mm;
}

// per-read QC and cap MQ
void compute_read_qc(const bam1_t *b, int C, ReadQC *out){
    memset(out,0,sizeof(*out)); out->mapq_raw=b->core.qual;
    const uint32_t *cig=bam_get_cigar(b); int nc=b->core.n_cigar;
    uint8_t *qual=bam_get_qual(b);
    int scL=left_clip_len(b,BAM_CSOFT_CLIP), scR=right_clip_len(b,BAM_CSOFT_CLIP);
    int hcL=left_clip_len(b,BAM_CHARD_CLIP), hcR=right_clip_len(b,BAM_CHARD_CLIP);
    int soft=scL+scR, hard=hcL+hcR;
    out->clipped_bases=soft+hard;

    long long sc_qsum=0;
    if (qual && scL){ for(int i=0;i<scL;i++) sc_qsum+=qual[i]; }
    if (qual && scR){ int lq=b->core.l_qseq; for(int i=lq-scR;i<lq;i++) sc_qsum+=qual[i]; }
    out->ClipQ = (int)(sc_qsum + 13LL*hard);

    int rpos=0;
    for(int k=0;k<nc;k++){
        int op=bam_cigar_op(cig[k]), ln=bam_cigar_oplen(cig[k]);
        if (op==BAM_CMATCH||op==BAM_CEQUAL||op==BAM_CDIFF){
            for(int i=0;i<ln;i++){ int q=qual?qual[rpos+i]:0; if (q>=BQ_MIN) out->M_q13++; }
            rpos+=ln;
        } else if (op==BAM_CINS){ out->ins_len+=ln; rpos+=ln; }
        else if (op==BAM_CDEL){ out->del_len+=ln; }
        else if (op==BAM_CSOFT_CLIP){ rpos+=ln; }
    }
    int mmn=0; Mismatch *mm = mismatches_from_MD(b,&mmn);
    for(int i=0;i<mmn;i++){
        if (mm[i].bq >= BQ_MIN){
            out->X_ge13++;
            out->SubQ += (mm[i].bq < BQ_CAP ? mm[i].bq : BQ_CAP);
        }
    }
    if (mm) free(mm);

    double t=0.0;
    if (out->X_ge13>0 && out->M_q13>0){
        double log10_t = out->X_ge13 * log10((double)out->M_q13)
                       - ( log10(tgamma((double)out->X_ge13 + 1.0)) / log(10.0) );
        t = (double)out->SubQ - 10.0*log10_t + (double)out->ClipQ/5.0;
    } else {
        t = (double)out->ClipQ/5.0;
    }
    if (t > C) out->cap_mq = -1;
    else { if (t < 0) t = 0; double cap = sqrt((C - t)/C) * C; out->cap_mq = (int)(cap + 0.499); }

    int cap_c = out->cap_mq < 0 ? 0 : out->cap_mq;
    int cap60 = (int)((cap_c * 60.0) / (double)C + 0.5);
    if (cap60 > 60) { cap60 = 60; }
    else if (cap60 < 0) { cap60 = 0; }
    out->cap_mq60 = cap60;
    out->eff_mq60 = (out->mapq_raw < out->cap_mq60) ? out->mapq_raw : out->cap_mq60;
}

// apply read to interval & add local features
long long apply_read_to_interval(const bam1_t *b, int iv_tid, int iv_start, int iv_end, const ReadQC *rq, Site *sites){
    const uint32_t *cig=bam_get_cigar(b); int nc=b->core.n_cigar;
    if (b->core.tid != iv_tid) return 0;
    int64_t ref=b->core.pos; int rpos=0; long long aligned_bp=0;
    uint8_t *seq = bam_get_seq(b);
    bool rev = (b->core.flag & BAM_FREVERSE) != 0;
    for(int k=0;k<nc;k++){
        int op=bam_cigar_op(cig[k]), ln=bam_cigar_oplen(cig[k]);
        if (op==BAM_CMATCH||op==BAM_CEQUAL||op==BAM_CDIFF){
            for(int i=0;i<ln;i++){
                int64_t rp=ref+i;
                if (rp>=iv_start && rp<iv_end){
                    int idx=(int)(rp - iv_start);
                    sites[idx].depth++;
                    vpush(&sites[idx].mq, rq->mapq_raw);
                    vpush(&sites[idx].capmq, rq->cap_mq60);
                    vpush(&sites[idx].effmq, rq->eff_mq60);
                    vpush(&sites[idx].subq, rq->SubQ);
                    vpush(&sites[idx].clipq, rq->ClipQ);
                    sites[idx].clip_bases_sum += rq->clipped_bases;

                    int nt = bam_seqi(seq, rpos+i) & 0xF;
                    switch(nt){
                        case 1: if (rev) sites[idx].nA_rev++; else sites[idx].nA_fwd++; break;
                        case 2: if (rev) sites[idx].nC_rev++; else sites[idx].nC_fwd++; break;
                        case 4: if (rev) sites[idx].nG_rev++; else sites[idx].nG_fwd++; break;
                        case 8: if (rev) sites[idx].nT_rev++; else sites[idx].nT_fwd++; break;
                        default: if (rev) sites[idx].nN_rev++; else sites[idx].nN_fwd++; break;
                    }

                    // histograms
                    int mqbin = rq->mapq_raw/10; if(mqbin<0) mqbin=0; if(mqbin>6) mqbin=6; sites[idx].hist_mq[mqbin]++;
                    int efbin = rq->eff_mq60/10; if(efbin<0) efbin=0; if(efbin>6) efbin=6; sites[idx].hist_eff[efbin]++;
                    double clipfrac = (double)rq->clipped_bases / (double)b->core.l_qseq;
                    if (clipfrac < 0.0) { clipfrac = 0.0; }
                    else if (clipfrac > 1.0) { clipfrac = 1.0; }
                    int cfbin = (int)floor(clipfrac*10.0); if (cfbin<0) cfbin=0; if (cfbin>9) cfbin=9;
                    sites[idx].hist_clipfrac[cfbin]++;
                }
            }
            int64_t segL=ref, segR=ref+ln;
            int64_t ovL = segL>iv_start?segL:iv_start;
            int64_t ovR = segR<iv_end?segR:iv_end;
            if (ovR>ovL) aligned_bp += (ovR-ovL);
            ref+=ln; rpos+=ln;
        } else if (op==BAM_CDEL){
            for(int i=0;i<ln;i++){
                int64_t rp=ref+i;
                if (rp>=iv_start && rp<iv_end){
                    int idx=(int)(rp - iv_start);
                    sites[idx].depth++;
                }
            }
            if (ref>=iv_start && ref<iv_end){
                int idx=(int)(ref - iv_start);
                sites[idx].del_len_sum += ln;
            }
            ref+=ln;
        } else if (op==BAM_CINS){
            if (ref>=iv_start && ref<iv_end){
                int idx=(int)(ref - iv_start);
                sites[idx].ins_len_sum += ln;
            }
            rpos+=ln;
        } else if (op==BAM_CREF_SKIP){ ref+=ln; }
        else if (op==BAM_CSOFT_CLIP){ rpos+=ln; }
    }
    return aligned_bp;
}

void add_local_mismatches(const bam1_t *b, int iv_tid, int iv_start, int iv_end, Site *sites){
    if (b->core.tid != iv_tid) return;
    const uint32_t *cig=bam_get_cigar(b); int nc=b->core.n_cigar;
    int64_t ref=b->core.pos; int rpos=0, cap=0;
    for(int k=0;k<nc;k++){
        int op=bam_cigar_op(cig[k]), ln=bam_cigar_oplen(cig[k]);
        if (op==BAM_CMATCH||op==BAM_CEQUAL||op==BAM_CDIFF) cap+=ln;
        else if (op==BAM_CINS||op==BAM_CSOFT_CLIP) rpos+=ln;
        else if (op==BAM_CDEL||op==BAM_CREF_SKIP) ref+=ln;
    }
    int64_t *refpos = cap? (int64_t*)malloc(cap*sizeof(int64_t)) : NULL;
    if (cap && !refpos){ perror("malloc"); exit(1); }
    ref=b->core.pos; rpos=0; int idx=0;
    for(int k=0;k<nc;k++){
        int op=bam_cigar_op(cig[k]), ln=bam_cigar_oplen(cig[k]);
        if (op==BAM_CMATCH||op==BAM_CEQUAL||op==BAM_CDIFF){
            for(int i=0;i<ln;i++) refpos[idx++]=ref+i;
            ref+=ln; rpos+=ln;
        } else if (op==BAM_CINS){ rpos+=ln; }
        else if (op==BAM_CDEL||op==BAM_CREF_SKIP){ ref+=ln; }
        else if (op==BAM_CSOFT_CLIP){ rpos+=ln; }
    }
    uint8_t *mdp=bam_aux_get(b,"MD");
    if (mdp){
        const char *md=bam_aux2Z(mdp); int ai=0; const char *p=md;
        while (*p){
            if (isdigit((unsigned char)*p)){ int v=0; while(isdigit((unsigned char)*p)){ v=v*10+(*p-'0'); ++p; } ai+=v; }
            else if (*p=='^'){ ++p; while(*p && isalpha((unsigned char)*p)) ++p; }
            else if (isalpha((unsigned char)*p)){
                if (ai<cap){
                    int64_t rp=refpos[ai];
                    if (rp>=iv_start && rp<iv_end){ int sidx=(int)(rp - iv_start); sites[sidx].mismatches++; }
                    ai++;
                }
                ++p;
            } else ++p;
        }
    }
    if (refpos) free(refpos);
}
