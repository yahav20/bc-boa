# bc-boa — Bounded-Cost Bi-Objective Search (BCP-BOA*)

Overview
--------
This repository implements BCP-BOA*, a bounded-cost variant of BOA* for bi-objective shortest-path search (SoCS 2022). BCP-BOA* finds Pareto-optimal solution paths whose two costs lie within given bounds (b1,b2). Multiple ordering functions are supported: lex1, lex2, sellex, min, max, avg.

Build
-----
Requirements: gcc and make.

From the repository root:

  make bcp_boastar

This produces the `bcp_boastar` executable (see `make all` for other targets).

Running bcp_boastar
-------------------
Usage:

  bcp_boastar [graph_file] [start_node] [goal_node] [ordering] [pivot] [zone]

Parameters:
- ordering: lex1 | lex2 | sellex | min | max | avg
- pivot   : ftl | tl | md | br | fbr
- zone    : 2 | 3 | 4

Example (BAY map, 1-based node indices):

  ./bcp_boastar Maps/BAY-road-d.txt 1 200 lex1 ftl 3

Notes for experiments
---------------------
- To compare ordering functions, run the same graph/start/goal/pivot/zone with different `ordering` values and compare the printed statistics (states generated/expanded, runtime).
- `sellex` automatically selects between lex1 and lex2 based on bounds.
- Maps are in `Maps/` (e.g., `BAY-road-d.txt`).
- For reproducible results, compile with optimizations (default CFLAGS include -O3) and redirect output to log files.

Key source files
----------------
- src/bcp_boastar.c — main BCP-BOA* implementation
- src/main_bc_boastar.c — CLI entry-point for running experiments
- src/node.h — graph node/search-node data structures
- Bounded-Cost_Bi-Objective_Search.md — paper and algorithm description

Example batch script (bash)
---------------------------
Run all orderings and save logs:

  for ord in lex1 lex2 sellex min max avg; do
    ./bcp_boastar Maps/BAY-road-d.txt 1 200 $ord ftl 3 | tee out_$ord.log
  done

Support
-------
If build or runtime errors occur, run `make clean` and rebuild. Paste compiler/linker output here and I’ll fix issues.
