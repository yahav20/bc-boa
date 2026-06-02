/////////////////////////////////////////////////////////////////////
// Carlos Hernandez & Yahav Ezer
// All rights reserved
/////////////////////////////////////////////////////////////////////

#include "include.h"
#include "bcp_boastar.h"
#include "graph.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static OrderingFunction parse_ordering(const char *s) {
    if (strcmp(s, "lex1")   == 0) return ORDER_LEX1;
    if (strcmp(s, "lex2")   == 0) return ORDER_LEX2;
    if (strcmp(s, "sellex") == 0) return ORDER_SEL_LEX;
    if (strcmp(s, "min")    == 0) return ORDER_MIN;
    if (strcmp(s, "max")    == 0) return ORDER_MAX;
    if (strcmp(s, "avg")    == 0) return ORDER_AVG;
    fprintf(stderr, "Unknown ordering '%s'. Use: lex1 lex2 sellex min max avg\n", s);
    exit(1);
}

static PivotType parse_pivot(const char *s) {
    if (strcmp(s, "ftl") == 0) return PIVOT_FTL;
    if (strcmp(s, "tl")  == 0) return PIVOT_TL;
    if (strcmp(s, "md")  == 0) return PIVOT_MD;
    if (strcmp(s, "br")  == 0) return PIVOT_BR;
    if (strcmp(s, "fbr") == 0) return PIVOT_FBR;
    fprintf(stderr, "Unknown pivot '%s'. Use: ftl tl md br fbr\n", s);
    exit(1);
}

/*----------------------------------------------------------------------------------*/
int main(int argc, char** argv) {
    if (argc != 7) {
        printf("Usage: %s [graph_file] [start_node] [goal_node] [ordering] [pivot] [zone]\n", argv[0]);
        printf("  ordering : lex1 | lex2 | sellex | min | max | avg\n");
        printf("  pivot    : ftl | tl | md | br | fbr\n");
        printf("  zone     : 2 | 3 | 4\n");
        exit(1);
    }

    char filename[128];
    strcpy(filename, argv[1]);
    start = atoi(argv[2]) - 1;
    goal  = atoi(argv[3]) - 1;

    OrderingFunction ord  = parse_ordering(argv[4]);
    PivotType        piv  = parse_pivot(argv[5]);
    int              zone = atoi(argv[6]);

    if (zone < 2 || zone > 4) {
        fprintf(stderr, "Zone must be 2, 3, or 4.\n");
        exit(1);
    }

    read_adjacent_table(filename);
    new_graph();
    call_bc_boastar(ord, piv, zone);
    return 0;
}
