#include <math.h>
#include "common.h"

static inline int clampi(int x, int a, int b){ return x<a?a:(x>b?b:x); }
static inline double clampd(double x, double a, double b){ return x<a?a:(x>b?b:x); }

void hist_accum_mq(int *hist7, int mq0_60){
    int b = mq0_60 / 10; if (b<0) b=0; if (b>6) b=6; hist7[b]++;
}
void hist_accum_eff(int *hist7, int eff0_60){
    int b = eff0_60 / 10; if (b<0) b=0; if (b>6) b=6; hist7[b]++;
}
void hist_accum_clipfrac(int *hist10, double f){
    int b = (int)floor(f*10.0); if (b<0) b=0; if (b>9) b=9; hist10[b]++;
}

void hist_ks_w1_js(const long long *hA, const long long *hB, int nbins, double binw,
                   double *ks, double *w1, double *js)
{
    long long NA=0, NB=0;
    for(int b=0;b<nbins;b++){ NA+=hA[b]; NB+=hB[b]; }
    *ks=NAN; *w1=NAN; *js=NAN;
    if (NA==0 || NB==0) return;
    double cA=0.0, cB=0.0, mjs=0.0, mks=0.0, mw1=0.0;
    for(int b=0;b<nbins;b++){
        double pA=(double)hA[b]/(double)NA, pB=(double)hB[b]/(double)NB;
        cA+=pA; cB+=pB;
        double d=fabs(cA-cB);
        if (d>mks) mks=d;
        mw1 += d * binw;
        double m=0.5*(pA+pB);
        if (pA>0 && m>0) mjs += 0.5*pA*log(pA/m);
        if (pB>0 && m>0) mjs += 0.5*pB*log(pB/m);
    }
    *ks = mks;
    *w1 = mw1;
    *js = sqrt(mjs / log(2.0));
}
