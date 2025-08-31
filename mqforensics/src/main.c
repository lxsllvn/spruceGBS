#include <stdio.h>
#include <string.h>
#include "common.h"

int analyze_main(int argc, char **argv);
int summarize_main(int argc, char **argv);

static void usage(void){
    fprintf(stderr,
    "mq_forensics <analyze|summarize> [options]\n"
    "  analyze   -b BAM -r BED -C INT -d MIN -o per_site.tsv -O per_interval.tsv [-f ref.fa] [--flank N] [--ref-only]  [--emit-suffstats] [--emit-hist]\n"
    "  summarize -i suffstats.sorted.tsv -o pooled.tsv\n");
}

int main(int argc, char **argv){
    if (argc < 2){ usage(); return 2; }
    if (!strcmp(argv[1], "analyze")) return analyze_main(argc-1, argv+1);
    if (!strcmp(argv[1], "summarize")) return summarize_main(argc-1, argv+1);
    usage();
    return 2;
}
