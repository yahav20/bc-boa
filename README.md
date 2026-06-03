# Bounded-Cost Bi-Objective Search: BCP-BOA* with Budget Score Ordering

A faithful implementation of **BCP-BOA*** (Skyler et al., SoCS 2022), extended with a novel
ordering function — **Budget Score (BS)** — proposed as part of an ongoing research project.

---

## Research Context

Multi-objective shortest-path search seeks solutions that are Pareto-optimal across two or more
cost dimensions simultaneously. **BOA*** (Bi-Objective A\*) solves this exactly but may expand
an exponential number of nodes when the Pareto-Optimal Frontier (POF) is large.

**BCP-BOA*** addresses this by operating in *bounded-cost* mode: given absolute cost bounds
(b₁, b₂), it finds a single solution satisfying both constraints, typically near a user-specified
*pivot point* on the POF. The key algorithmic choices — how to order the open list and how to
prune dominated nodes — dramatically affect runtime and solution quality.

This project studies the effect of different ordering functions on BCP-BOA* performance,
and introduces a new ordering function designed to exploit the budget structure of the problem.

---

## Novel Contribution: Budget Score (BS) Ordering

### Motivation

Existing ordering functions (Lex, Min, Max, Avg) either use unnormalised costs, which bias
search toward one objective, or use range-normalised costs, which are query-independent and
ignore the actual budget constraints (b₁, b₂).

### Definition

Let the *budget-normalised f-values* be:

```
F₁(n) = f₁(n) / b₁      F₂(n) = f₂(n) / b₂
```

Both F₁ and F₂ lie in [0, 1] for any node that survives bound pruning.
The dynamic weight is:

```
aₙ = F₁(n) / (F₁(n) + F₂(n))
```

which is large when the node is relatively more constrained in dimension 1, and small otherwise.
The **Budget Score** is:

```
S(n) = aₙ · F₁(n) + (1 − aₙ) · F₂(n)  =  (F₁² + F₂²) / (F₁ + F₂)
```

This is the weighted sum of the two normalised f-values, where each objective is weighted by its
own relative pressure. Geometrically, S(n) is minimised along the diagonal of the feasible
budget square, naturally guiding the search toward *balanced* solutions.

### Dominance Strategy: Relaxed Single-Label Pruning

Classical dominance (used by Min/Max/Avg) maintains a Pareto front of all expanded (g₁, g₂)
pairs per node — O(|POF|) space and O(|POF|) time per check.

For BS we propose **relaxed scalar dominance**: store only the minimum S value seen so far at
each graph node:

```
g_best_S[v] = min { S(n) : n expanded at state v, current round }
```

A new node at state v with score S ≥ g_best_S[v] is pruned. This reduces per-node dominance
from O(|POF|) to **O(1)** in both time and space, and enables early pruning at *generation*
time (not just expansion time).

**Trade-off**: relaxed dominance is not Pareto-complete — it may discard a (g₁, g₂) pair that
is non-dominated in the traditional sense but has a higher S score than a previously expanded
path. The algorithm is still sound (every returned solution satisfies the bounds), but it operates
as a *greedy* approximation to bounded-cost Pareto search.

### Solution Quality Metric

Because BS may return a solution that is not strictly Pareto-optimal, the benchmark reports a
quality metric for every BS result. Given the BS solution (c₁, c₂) and the full POF computed
by BOA*, the nearest POF point (p₁, p₂) is found by normalised Euclidean distance:

```
dist = √( ((c₁−p₁)/(max₁−min₁))² + ((c₂−p₂)/(max₂−min₂))² )
```

The relative error is then:

```
Error = (c₁ − p₁)/p₁ + (c₂ − p₂)/p₂
```

An error of 0.00 means BS found an exact Pareto-optimal solution.

---

## Ordering Functions

| Name    | Key formula | Dominance | Complexity |
|---------|-------------|-----------|------------|
| Lex1    | f₁·BASE + f₂ | scalar g₂min | O(1) |
| Lex2    | f₂·BASE + f₁ | scalar g₁min | O(1) |
| Sel-Lex | Lex1 or Lex2 (by b̄₁ vs b̄₂) | scalar | O(1) |
| Min     | min(nf₁,nf₂)·BASE + max | Pareto front | O(\|POF\|) |
| Max     | max(nf₁,nf₂)·BASE + min | Pareto front | O(\|POF\|) |
| Average | avg(nf₁,nf₂)·BASE + min | Pareto front | O(\|POF\|) |
| **BS**  | **(F₁²+F₂²)/(F₁+F₂)·BASE + min(F₁,F₂)** | **scalar S** | **O(1)** |

nf₁, nf₂ = range-normalised f-values; F₁, F₂ = budget-normalised f-values.

---

## Benchmark Protocol

The benchmark replicates the experimental setup of Skyler et al., SoCS 2022.

For each query the benchmark runs:

1. **Heuristic computation** — backward Dijkstra on both cost dimensions.
2. **Map-extreme computation** — forward Lex Dijkstra to find (min₁, max₁, min₂, max₂).
3. **Full BOA*** — to obtain the complete Pareto-Optimal Frontier and select pivot points.
4. **BCP-BOA*** — for every combination of:
   - 7 ordering functions: Lex1, Lex2, Sel-Lex, Min, Max, Average, **BS**
   - 5 pivot points: FTL, TL, MD, BR, FBR
   - 4 zones: Z2, Z3, Z4, Z5
   - Total: **112 combinations per query**

Summary tables report average expansions, search time, states pruned, states generated, and
success rate. For BS an additional quality table reports the average solution error.

---

## Build

```bash
make            # builds all targets: boa, bcp_boastar, benchmark
make clean      # remove binaries and object files
```

Requirements: `gcc`, `make`.

---

## Usage

### Single query

```bash
./bcp_boastar <map_file> <start> <goal> <ordering> <pivot> <zone>

# ordering: lex1 | lex2 | sellex | min | max | avg | bs
# pivot   : ftl | tl | md | br | fbr
# zone    : 2 | 3 | 4

./bcp_boastar Maps/BAY-road-d.txt 217950 116998 bs md 3
```

### Full benchmark

```bash
./benchmark <map_file> <num_queries> <timeout_sec>

./benchmark Maps/BAY-road-d.txt 50 30
```

Maps and query files must reside in `Maps/` and `Queries/` respectively.
The benchmark derives the query filename automatically from the map filename prefix
(e.g. `BAY-road-d.txt` → `Queries/BAY-queries`).

---

## Repository Structure

```
src/
  bcp_boastar.c / .h   — BCP-BOA* core: ordering, dominance, budget score
  main_bc_boastar.c    — single-query CLI
  benchmark.c          — full benchmark replicating SoCS 2022 setup
  boastar.c / .h       — exact BOA* (used for POF computation and heuristics)
  node.h               — graph node and search node data structures
  heap.c / .h          — binary heap (open list)
  graph.c / .h         — graph loading and adjacency tables
  pool.c               — malloc-based node allocator (single-query binary)
  pool_bench.c         — slab-based node allocator (benchmark binary)
Maps/                  — road network graphs (BAY, NY, …)
Queries/               — query files (start/goal pairs)
Results/               — experiment output logs
```

---

## Reference

> T. Skyler, C. Hernandez, J. Baier, S. Koenig.
> *Bounded-Cost Bi-Objective Search.*
> Proceedings of the 15th Annual Symposium on Combinatorial Search (SoCS), 2022.
