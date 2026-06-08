# Results — Benchmark Scripts and Output

This directory contains the benchmark runner, aggregation tool, and experiment output
for the BCP-BOA* evaluation.

---

## Directory Contents

| File | Description |
|------|-------------|
| `run_chunks.sh` | Shell script: runs the benchmark in sequential chunks, saves log and CSV |
| `aggregate.py` | Python script: combines per-chunk CSV data into a single summary table |
| `BAY_chunks.csv` | Raw per-chunk accumulator data for the BAY road network |
| `NY_chunks.csv` | Raw per-chunk accumulator data for the NY road network |
| `bay_results.txt` | Aggregated summary tables for BAY (49 queries) |
| `ny_results.txt` | Aggregated summary tables for NY |

---

## Running the Benchmark

### Option A — Direct run (few queries)

```bash
# from the repository root
./benchmark <map_file> <num_queries> <timeout_sec> [start_query]
```

Example:
```bash
./benchmark Maps/BAY-road-d.txt 10 30
```

Runs 10 queries, 30 s timeout per ordering, prints a full summary table, and
appends raw data to `BAY_chunks.csv`.  The optional `start_query` argument (default 0)
sets the 0-based offset into the query file, so successive partial runs can cover
disjoint ranges.

### Option B — Chunked run (recommended for many queries)

`run_chunks.sh` splits a large run into smaller sequential chunks.  Each chunk calls
`./benchmark` with the appropriate `start_query` offset, so the pool is freed between
chunks and intermediate results are visible as the run progresses.

```bash
# from the repository root
./Results/run_chunks.sh <map_file> <total_queries> <timeout_sec> [chunk_size=5]
```

| Argument | Description | Default |
|----------|-------------|---------|
| `map_file` | Path to the road network graph | — |
| `total_queries` | Total number of queries to run | — |
| `timeout_sec` | Per-ordering time limit in seconds | — |
| `chunk_size` | Queries per chunk | `5` |

Example — run all 50 BAY queries in chunks of 5, 30 s timeout:
```bash
./Results/run_chunks.sh Maps/BAY-road-d.txt 50 30 5
```

**Outputs:**
- Console output is mirrored to `Results/BAY_all.log`
- Raw cumulative sums are appended to `BAY_chunks.csv` after each chunk

### Option C — Aggregate after a chunked run

```bash
python3 Results/aggregate.py BAY_chunks.csv
```

Sums all chunk rows and prints the same summary tables that `./benchmark` prints,
averaged over the total number of queries in the CSV.

---

## Output Tables

Both `./benchmark` and `aggregate.py` print the following tables.

### BOA* Baseline

Average cost of the shared setup phase (run once per query, before any BCP-BOA* call):

| Metric | Meaning |
|--------|---------|
| Setup time (ms) | Backward Dijkstra ×2 + Lex Dijkstra ×2 + full BOA* |
| States expanded | BOA* node expansions |
| States generated | Successors added to the open list |
| POF solutions | Number of Pareto-optimal solutions found |

### BCP-BOA* Tables

Five tables — **Expansions**, **Search Time**, **Pruned**, **Success Rate**,
**Generated** — each with rows for every ordering function and columns for every
(zone, pivot) combination.

#### Columns: Zone × Pivot

Zones set how tight the budget bounds are relative to the pivot point.
All bounds are expressed as normalised fractions of the cost range [min_c, max_c].

| Zone | b̄₁, b̄₂ formula (given pivot p) | Tightness |
|------|----------------------------------|-----------|
| Z2 | (3p+1)/4 | Tightest — very close to the pivot |
| Z3 | (p+1)/2 | Medium |
| Z4 | (p+3)/4 | Loose |
| Z5 | 1.0, 1.0 | Loosest — full cost range |

Pivots select where on the Pareto-Optimal Frontier (POF) the budget is centred.
Coordinates are normalised to [0,1]² by the range of optimal costs.

| Pivot | Location |
|-------|----------|
| FTL | Far Top-Left: Lex-optimal in dimension 1 (min c₁, coordinates (0, 1)) |
| TL | Top-Left: closest POF point to midpoint(FTL, MD) |
| MD | Middle: closest POF point to the normalised diagonal |
| BR | Bottom-Right: closest POF point to midpoint(FBR, MD) |
| FBR | Far Bottom-Right: Lex-optimal in dimension 2 (min c₂, coordinates (1, 0)) |

#### Rows: Ordering Functions

| Ordering | Key formula | Dominance | Complexity |
|----------|-------------|-----------|------------|
| Lex1 | f₁·BASE + f₂ | scalar g₂min | O(1) |
| Lex2 | f₂·BASE + f₁ | scalar g₁min | O(1) |
| Sel-Lex | Lex1 if b̄₁ ≤ b̄₂, else Lex2 | scalar | O(1) |
| Min | min(nf₁,nf₂)·BASE + max | Pareto front | O(\|POF\|) |
| Max | max(nf₁,nf₂)·BASE + min | Pareto front | O(\|POF\|) |
| Average | avg(nf₁,nf₂)·BASE + min | Pareto front | O(\|POF\|) |
| **BS** | **(F₁²+F₂²)/(F₁+F₂)·BASE + min(F₁,F₂)** | scalar S-score | O(1) |

nf = range-normalised f-value; F = budget-normalised f-value (F = f/b).

### BS Solution Quality

Because BS uses relaxed scalar dominance it may return a solution that is not
strictly Pareto-optimal.  This table reports the average **relative error** versus
the nearest POF point found by BOA*:

```
Error = (c₁ − p₁)/p₁ + (c₂ − p₂)/p₂
```

where (c₁, c₂) is the BS solution and (p₁, p₂) is the nearest POF point by
normalised Euclidean distance.

- `0.0000` — BS found an exact Pareto-optimal solution.
- `N/A` — no solution found (timeout or bounds infeasible).

---

## Raw CSV Format

`<MAP>_chunks.csv` stores raw sums and counts so that multiple runs can be combined
without precision loss.  Two row types are written by `./benchmark`:

```
# BOA* baseline row (one per chunk):
boa_stats,<map>,<start_q>,<n_queries>,<boa_exp_sum>,<boa_gen_sum>,<boa_sol_sum>,<setup_time_sum>

# BCP-BOA* row (one per (ordering, zone, pivot) combination per chunk):
bcp_stats,<map>,<start_q>,<ord_i>,<zone_i>,<piv_i>,<exp_sum>,<gen_sum>,<prune_sum>,
          <time_sum>,<bbar1_sum>,<bbar2_sum>,<solved>,<timed_out>,<n>,<error_sum>,<error_n>
```

Index encodings: `ord_i` ∈ {0=Lex1, 1=Lex2, 2=Sel-Lex, 3=Min, 4=Max, 5=Average, 6=BS},
`zone_i` ∈ {0=Z2, 1=Z3, 2=Z4, 3=Z5}, `piv_i` ∈ {0=FTL, 1=TL, 2=MD, 3=BR, 4=FBR}.
