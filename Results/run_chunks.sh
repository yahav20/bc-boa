#!/usr/bin/env bash
# run_chunks.sh — Run benchmark in small chunks and save output.
#
# Usage:
#   ./run_chunks.sh <map_file> <total_queries> <timeout_sec> [chunk_size=5]
#
# Example:
#   ./run_chunks.sh Maps/BAY-road-d.txt 20 30 5
#
# Each chunk runs `chunk_size` consecutive queries starting where the
# previous chunk left off. Output from all chunks is appended to
#   results/<MAP>_all.log
# Raw sums are accumulated in
#   <MAP>_chunks.csv   (written by the benchmark binary itself)
# After all chunks finish, run:
#   python3 aggregate.py <MAP>_chunks.csv

set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 <map_file> <total_queries> <timeout_sec> [chunk_size=5]"
    exit 1
fi

MAP="$1"
TOTAL="$2"
TIMEOUT="$3"
CHUNK="${4:-5}"

if [[ ! -f "$MAP" ]]; then
    echo "Error: map file '$MAP' not found."
    exit 1
fi

if [[ ! -x "./benchmark" ]]; then
    echo "Error: ./benchmark not found or not executable. Run 'make benchmark' first."
    exit 1
fi

# Derive the map prefix (e.g. "BAY" from "Maps/BAY-road-d.txt")
BASENAME=$(basename "$MAP")
PREFIX="${BASENAME%%-*}"

OUTDIR="results"
mkdir -p "$OUTDIR"
LOG="$OUTDIR/${PREFIX}_all.log"
CSV="${PREFIX}_chunks.csv"

# Fresh run: clear previous log and CSV
rm -f "$LOG" "$CSV"

echo "=========================================="
echo " BCP-BOA* Chunked Benchmark"
echo " Map     : $MAP"
echo " Queries : $TOTAL  |  Chunk size: $CHUNK"
echo " Timeout : ${TIMEOUT}s per ordering"
echo " Log     : $LOG"
echo " CSV     : $CSV"
echo "=========================================="
echo ""

start=0
chunk_num=1
while (( start < TOTAL )); do
    remaining=$(( TOTAL - start ))
    count=$(( remaining < CHUNK ? remaining : CHUNK ))
    end=$(( start + count ))

    echo "=========================================="  | tee -a "$LOG"
    echo " Chunk $chunk_num: queries $((start+1))-${end} of $TOTAL" | tee -a "$LOG"
    echo "==========================================" | tee -a "$LOG"

    ./benchmark "$MAP" "$count" "$TIMEOUT" "$start" 2>&1 | tee -a "$LOG"

    echo "" | tee -a "$LOG"

    start=$(( start + count ))
    chunk_num=$(( chunk_num + 1 ))
done

echo "=========================================="
echo " All $chunk_num chunks complete."
echo " Full log : $LOG"
echo " Raw CSV  : $CSV"
echo ""
echo " To produce a combined summary table run:"
echo "   python3 aggregate.py $CSV"
echo "=========================================="
