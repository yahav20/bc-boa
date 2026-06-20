#!/usr/bin/env python3
"""
aggregate.py — Combine per-chunk CSV files into a single summary table.

Usage:
    python3 aggregate.py <map>_chunks.csv

The CSV is written by ./benchmark (with a start_query argument).
Each run appends its raw sums+counts so multiple chunks accumulate in one file.
This script sums all chunks and prints the same tables that benchmark.c prints.
"""

import sys
import csv
from collections import defaultdict

ORD_NAMES  = ["Lex1", "Lex2", "Sel-Lex", "Min", "Max", "Average", "BS", "SBS"]
PIV_NAMES  = ["FTL", "TL", "MD", "BR", "FBR"]
ZONES      = [2, 3, 4, 5]
N_ORD      = 8
N_ZONES    = 4
N_PIV      = 5

# ---- Data containers ----

class BoaAcc:
    def __init__(self):
        self.n              = 0
        self.boa_exp_sum    = 0.0
        self.boa_gen_sum    = 0.0
        self.boa_sol_sum    = 0.0
        self.setup_time_sum = 0.0

class BcpAcc:
    def __init__(self):
        self.n          = 0
        self.exp_sum    = 0.0
        self.gen_sum    = 0.0
        self.prune_sum  = 0.0
        self.time_sum   = 0.0
        self.bbar1_sum  = 0.0
        self.bbar2_sum  = 0.0
        self.solved     = 0
        self.timed_out  = 0
        self.error_sum  = 0.0
        self.error_n    = 0

# ---- Load CSV ----

def load(path):
    boa = BoaAcc()
    bcp = defaultdict(BcpAcc)   # key: (ord_i, zone_i, piv_i)

    with open(path, newline="") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(",")
            rtype = parts[0]

            if rtype == "boa_stats":
                # type,map,start_q,n,boa_exp_sum,boa_gen_sum,boa_sol_sum,setup_time_sum
                boa.n              += int(parts[3])
                boa.boa_exp_sum    += float(parts[4])
                boa.boa_gen_sum    += float(parts[5])
                boa.boa_sol_sum    += float(parts[6])
                boa.setup_time_sum += float(parts[7])

            elif rtype == "bcp_stats":
                # type,map,start_q,ord_i,zone_i,piv_i,exp_sum,gen_sum,prune_sum,
                #   time_sum,bbar1_sum,bbar2_sum,solved,timed_out,n,error_sum,error_n
                key = (int(parts[3]), int(parts[4]), int(parts[5]))
                a = bcp[key]
                a.exp_sum   += float(parts[6])
                a.gen_sum   += float(parts[7])
                a.prune_sum += float(parts[8])
                a.time_sum  += float(parts[9])
                a.bbar1_sum += float(parts[10])
                a.bbar2_sum += float(parts[11])
                a.solved    += int(parts[12])
                a.timed_out += int(parts[13])
                a.n         += int(parts[14])
                a.error_sum += float(parts[15])
                a.error_n   += int(parts[16])

    return boa, bcp

# ---- Printing helpers ----

COL_W = 8
LBL_W = 10

def hline(n):
    print("-" * n)

def print_table(title, metric, bcp):
    ncols = N_PIV * 3 + 1   # 16 columns
    print(f"\n--- {title} ---")

    # Header row
    row = f"{'Ordering':<{LBL_W}}"
    for zi in range(3):
        for pi in range(N_PIV):
            label = f"Z{ZONES[zi]}:{PIV_NAMES[pi]}"
            row += f" {label:>{COL_W}}"
    row += f" {'Z5:Any':>{COL_W}}"
    print(row)

    # b1 avg row
    row = f"{'b1(avg)':<{LBL_W}}"
    for zi in range(3):
        for pi in range(N_PIV):
            a = bcp.get((0, zi, pi))
            avg = a.bbar1_sum / a.n if a and a.n > 0 else 0.0
            row += f" {avg:{COL_W}.2f}"
    a = bcp.get((0, 3, 0))
    avg = a.bbar1_sum / a.n if a and a.n > 0 else 1.0
    row += f" {avg:{COL_W}.2f}"
    print(row)

    # b2 avg row
    row = f"{'b2(avg)':<{LBL_W}}"
    for zi in range(3):
        for pi in range(N_PIV):
            a = bcp.get((0, zi, pi))
            avg = a.bbar2_sum / a.n if a and a.n > 0 else 0.0
            row += f" {avg:{COL_W}.2f}"
    a = bcp.get((0, 3, 0))
    avg = a.bbar2_sum / a.n if a and a.n > 0 else 1.0
    row += f" {avg:{COL_W}.2f}"
    print(row)

    hline(LBL_W + ncols * (COL_W + 1))

    for oi in range(N_ORD):
        row = f"{ORD_NAMES[oi]:<{LBL_W}}"
        for zi in range(3):
            for pi in range(N_PIV):
                v = metric.get((oi, zi, pi))
                if v is not None:
                    row += f" {v:{COL_W}.1f}"
                else:
                    row += f" {'N/A':>{COL_W}}"
        v = metric.get((oi, 3, 0))
        if v is not None:
            row += f" {v:{COL_W}.1f}"
        else:
            row += f" {'N/A':>{COL_W}}"
        print(row)

def build_metrics(bcp):
    avg_exp   = {}
    avg_time  = {}
    avg_prune = {}
    avg_gen   = {}
    pct_ok    = {}

    for (oi, zi, pi), a in bcp.items():
        if a.n > 0:
            avg_exp  [(oi,zi,pi)] = a.exp_sum   / a.n / 1000.0
            avg_time [(oi,zi,pi)] = a.time_sum  / a.n
            avg_prune[(oi,zi,pi)] = a.prune_sum / a.n / 1000.0
            avg_gen  [(oi,zi,pi)] = a.gen_sum   / a.n / 1000.0
            pct_ok   [(oi,zi,pi)] = 100.0 * a.solved / a.n

    return avg_exp, avg_time, avg_prune, avg_gen, pct_ok

# ---- Main ----

def main():
    if len(sys.argv) < 2:
        print(f"Usage: python3 {sys.argv[0]} <map>_chunks.csv")
        sys.exit(1)

    path = sys.argv[1]
    boa, bcp = load(path)

    if boa.n == 0:
        print("No data found in CSV.")
        sys.exit(1)

    N = float(boa.n)
    print("=================================================================")
    print(f"AGGREGATED SUMMARY ({boa.n} queries total, from {path})")
    print("=================================================================")

    print(f"\n--- BOA* BASELINE (avg over {boa.n} queries) ---")
    print(f"  Setup (heuristic+BOA*) time : {boa.setup_time_sum/N:7.2f} ms")
    print(f"  BOA* states expanded        : {boa.boa_exp_sum/N:7.0f}")
    print(f"  BOA* states generated       : {boa.boa_gen_sum/N:7.0f}")
    print(f"  POF solutions found         : {boa.boa_sol_sum/N:7.1f}")

    avg_exp, avg_time, avg_prune, avg_gen, pct_ok = build_metrics(bcp)

    print_table("EXPANSIONS (avg, thousands)", avg_exp, bcp)
    print_table("SEARCH TIME (avg, ms)",       avg_time, bcp)
    print_table("PRUNED (avg, thousands)",      avg_prune, bcp)
    print_table("SUCCESS RATE (%)",             pct_ok, bcp)
    print_table("GENERATED (avg, thousands)",   avg_gen, bcp)

    # BS/SBS solution quality
    ORDER_BS  = 6
    ORDER_SBS = 7
    print("\n--- BS/SBS SOLUTION QUALITY: avg error = (c1-p1)/p1 + (c2-p2)/p2 ---")
    print("  (p1,p2) = nearest POF point by normalised Euclidean distance)")
    print("  0.00 = exact Pareto point found;  N/A = no solution\n")

    row = f"{'Zone':<8}"
    for zi in range(3):
        for pi in range(N_PIV):
            label = f"Z{ZONES[zi]}:{PIV_NAMES[pi]}"
            row += f" {label:>7}"
    row += f" {'Z5:Any':>8}"
    print(row)
    hline(8 + (3*N_PIV+1)*8)

    for oi, label in ((ORDER_BS, "BS"), (ORDER_SBS, "SBS")):
        row = f"{label:<8}"
        for zi in range(3):
            for pi in range(N_PIV):
                a = bcp.get((oi, zi, pi))
                if a and a.error_n > 0:
                    row += f" {a.error_sum/a.error_n:7.4f}"
                else:
                    row += f" {'N/A':>7}"
        a = bcp.get((oi, 3, 0))
        if a and a.error_n > 0:
            row += f" {a.error_sum/a.error_n:8.4f}"
        else:
            row += f" {'N/A':>8}"
        print(row)
    print()

if __name__ == "__main__":
    main()
