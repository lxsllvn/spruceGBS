#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "common.h"

void vpush(vecd *x, double a){
    if (x->n == x->cap){
        x->cap = x->cap ? x->cap*2 : 32;
        x->v = (double*)realloc(x->v, x->cap * sizeof(double));
        if (!x->v){ perror("realloc"); exit(1); }
    }
    x->v[x->n++] = a;
}
static int dcmp(const void *a, const void *b){
    double x = *(const double*)a, y = *(const double*)b;
    return (x<y)?-1:(x>y);
}
double dmean(const vecd *x){
    if (!x->n) return 0.0;
    double s=0; for(int i=0;i<x->n;i++) s+=x->v[i]; return s/x->n;
}
double dmedian(vecd *x){
    if (!x->n) return 0.0;
    qsort(x->v, x->n, sizeof(double), dcmp);
    if (x->n & 1) return x->v[x->n/2];
    return 0.5*(x->v[x->n/2-1] + x->v[x->n/2]);
}
double dsd(const vecd *x){
    if (x->n < 2) return 0.0;
    double mu=dmean(x), s2=0; for(int i=0;i<x->n;i++){ double d=x->v[i]-mu; s2+=d*d; }
    return sqrt(s2/(x->n-1));
}

// hist median with uniform within-bin interpolation
double median_from_hist_uniform_bins(const long long *h, int nbins, double bin0_left, double bin_width){
    long long tot=0; for(int b=0;b<nbins;b++) tot+=h[b];
    if (tot==0) return NAN;
    long long target=(tot+1)/2, cum=0;
    for(int b=0;b<nbins;b++){
        long long next=cum+h[b];
        if (target<=next){
            long long inbin = target - cum;
            double frac = (h[b]>0) ? ((double)inbin - 0.5) / (double)h[b] : 0.5;
            double left = bin0_left + b*bin_width;
            double right = left + bin_width;
            if (b==nbins-1) right = bin0_left + nbins*bin_width;
            return left + frac * (right-left);
        }
        cum=next;
    }
    return bin0_left + (nbins-0.5)*bin_width;
}
