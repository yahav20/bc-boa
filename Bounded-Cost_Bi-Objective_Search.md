# Bounded-Cost Bi-Objective Heuristic Search

**Proceedings of the Fifteenth International Symposium on Combinatorial Search (SoCS 2022)**

**Authors:**
* Shawn Skyler¹, Dor Atzmon¹, Ariel Felner¹, Oren Salzman²
* Han Zhang³, Sven Koenig³, William Yeoh⁴, Carlos Hernández Ulloa⁵

**Affiliations:**
1. Ben-Gurion University of the Negev
2. Technion - Israel Institute of Technology
3. University of Southern California
4. Washington University in St. Louis
5. Universidad Andrés Bello

**Contact:** shawn@post.bgu.ac.il, dorat@post.bgu.ac.il, felner@bgu.ac.il, osalzman@cs.technion.ac.il, zhan645@usc.edu, skoenig@usc.edu, wyeoh@wustl.edu, carlos.hernandez@uss.cl

---

## Abstract
There are many settings that extend the basic shortest-path search problem. In Bounded-Cost Search, we are given a constant bound, and the task is to find a solution within the bound. In Bi-Objective Search, each edge is associated with two costs (objectives), and the task is to minimize both objectives. In this paper, we combine both settings into a new setting of Bounded-Cost Bi-Objective Search. We are given two bounds, one for each objective, and the task is to find a solution within these bounds. We provide a scheme for normalizing the two objectives, introduce several algorithms for this new setting and compare them experimentally.

## Introduction and Background
A* (Hart, Nilsson, and Raphael 1968) and its many variants are used to solve classical shortest-path problems optimally. Nodes $n$ are expanded according to $f(n)=g(n)+h(n)$. If $h(n)$ is admissible, an optimal solution path will be found. Nevertheless, there are other settings that do not require optimal solution paths. One such setting is **Bounded-Cost Search (BCS)** (Stern et al. 2014). In BCS, we are given a bound $C$, and the task is to quickly find a solution path with cost $\le C$.

Another setting is **Bi-Objective Search (BOS)** (Raith and Ehrgott 2009). In BOS, we are given a graph (state-space) where two types of costs are associated with its edges that need to be minimized, e.g., travel distance and time for transportation problems. BOS has many applications, ranging from robotics (Fu et al. 2019; Fu, Salzman, and Alterovitz 2021) to transportation (Bronfman et al. 2015). The notion of an optimal path does not exist here. Instead, the "best" possible paths are those that are undominated by other paths. Path $\pi$ is said to dominate path $\pi^{\prime}$ iff both of $\pi$'s costs are not larger than the corresponding costs of $\pi^{\prime}$ and at least one cost of $\pi$ is strictly smaller than the corresponding cost of $\pi^{\prime}$. The **Pareto-Optimal Frontier (POF)** is the set of solution paths that are undominated by other solution paths. The BOA* algorithm (Hernández Ulloa et al. 2020) is considered the state-of-the-art algorithm for finding the POF.

In this paper we combine both BCS and BOS to define the problem of **Bounded-Cost Bi-Objective Search (BC-BOS)** where we are given a pair of bounds, one for each of the two costs, and need to find a solution path whose costs are below these bounds. In particular, in this paper, we are interested in the more restrictive problem of finding a solution path in the POF whose costs are below the bounds.

We first provide a normalization mechanism that converts the costs of both objectives (which could have different scales, e.g., time and distance) to values in the range of $[0...1]$. This places the two objectives on the same scale and, therefore, makes them directly comparable. We then introduce the BCP-BOA* algorithm, which modifies BOA* to return a POF solution path within the given bounds. BCP-BOA* uses our new notion of bound pruning, which prunes nodes that are not within the bounds. We introduce several variants of BCP-BOA* that explore the set of required nodes by using various ordering functions to order nodes in OPEN. We then experimentally compare the variants of BCP-BOA* on different pairs of cost bounds and study the advantages and disadvantages of each algorithm. Indeed, all variants run much faster than BOA*, that generates the full set of POF solution paths.

## Definitions and Background
A Bi-Objective Search (BOS) problem is characterized by a graph $G=(V,E)$, where $V$ is a set of vertices (states) and $E$ is a set of edges, a start vertex $start \in V$ and a goal vertex $goal \in V$. We use a boldface font to indicate a pair of two numbers, e.g., $\mathbf{b}=(b_1,b_2)$. The addition of two pairs $\mathbf{x}$ and $\mathbf{y}$ yields $(x_1+y_1, x_2+y_2)$. We say that $\mathbf{x}$ dominates $\mathbf{y}$ ($\mathbf{x} < \mathbf{y}$) iff $(x_1 \le y_1 \wedge x_2 < y_2)$ or $(x_1 < y_1 \wedge x_2 \le y_2)$. 

Each edge $e \in E$ is associated with two cost functions $c(e)=(c_1(e),c_2(e))$, i.e., a cost for each objective. A path $\pi=[v_1,...,v_n]$ is a list of neighboring vertices. The cost of path $\pi$ is $C(\pi)=(C_1(\pi),C_2(\pi))=\sum_{i=1}^{n-1}c(v_i,v_{i+1})$. A solution path is a path from $start$ to $goal$. Since there are two costs in BOS, the notion of an optimal solution path does not exist. Instead, the Pareto-Optimal Frontier (POF) is the set of solution paths $\Pi$ such that each path $\pi \in \Pi$ is not dominated by any other solution path.

For BOS, an admissible heuristic function $h(n)=(h_1(n),h_2(n))$ is a lower bound on the costs of a path from a given node $n$ to $goal$. We denote the heuristic for the node of the start vertex by $h(start)$ (resp. $h(goal)$ for the goal vertex). We also assume that $h$ is consistent. Specifically, we follow the common practice in BOS and use the individual shortest path heuristic function (Hernández Ulloa et al. 2020; Goldin and Salzman 2021; Pulido, Mandow, and Pérez-de-la Cruz 2015). That is, for $h(n)$, $h_i(n)$ is the cost-minimal path from $n$ to $goal$ using the $i$-th objective only.

## Previous Work on BOS
In BOS, the number of nodes can be exponentially larger than the number of vertices of the underlying graph because every path to a vertex has its unique node with possibly unique f-values. In fact, the POF itself can contain an exponential number of paths. Several A*-based search algorithms have been designed to find the POF. Examples include Multi-Objective A* (MOA*) (Stewart and White III 1991), NAMOA* (Mandow and De La Cruz 2005) and NAMOA*-dr (Pulido, Mandow, and Pérez-de-la Cruz 2015). All these algorithms use the same general best-first search framework, denoted here by BO-BFS. 

BOA* (Hernández Ulloa et al. 2020) is a BO-BFS algorithm, which is considered the current state-of-the-art. BOA* exploits the fact that nodes are ordered in lexicographic order of their f-values to perform all dominance checks in constant time. 

## Bounded-Cost BOS
The Bounded-Cost Bi-Objective Search (BC-BOS) problem is a BOS problem which is also characterized by a pair of bounds (budgets) $\mathbf{b}=(b_1,b_2)$. A solution path is now a path $\pi=[start,...,goal]$ whose costs are within the bounds, i.e., $C(\pi) \le \mathbf{b}$. The task is to find a solution path as fast as possible. In this paper, we focus on the more restrictive problem of finding a solution path that is also in the POF. We call this restricted problem BCP-BOS to differentiate it from the general BC-BOS problem.

## Normalizing the Two Objectives
In many problem instances, the two objectives may not be measured on the same scale, e.g., time in minutes and distance in miles. Therefore, we introduce a normalizing procedure that maps all values into the range $[0...1]$, which allows us to compare and combine both objectives since they are now on the same scale.

### Extreme Costs of Solution Paths
Let $min_1 = h_1(start)$ (the $c_1$-cost of the cost-minimal path from start to goal), and let $OPT_1$ be the set of all solution paths with cost $min_1$. Now, let $max_2 = \min_{\pi \in OPT_1} C_2(\pi)$ (namely, $max_2$ is the minimal $c_2$-cost of paths whose $c_1$-cost is $min_1$). $min_2$, $OPT_2$ and $max_1$ are defined analogously. 

$(min_1, max_2)$ and $(max_1, min_2)$ are the most extreme solution costs in the POF, and $h(start)=(min_1, min_2)$. Given $min_i$ and $max_i$, we can easily solve the BCP-BOS problem in the following cases: 
1. If one of the bounds satisfies $b_i < min_i$ then no solution path exists. 
2. If one of the bounds satisfies $b_i > max_i$ then one of the extreme solution paths in the POF is a valid solution. 
Thus, henceforth, we limit our discussion to the setting where $min_i \le b_i \le max_i$ for both $i \in \{1,2\}$.

### Normalizing the Costs
We now define a normalization function for any cost value $x_i$ with $min_i \le x_i \le max_i$ (for $i \in \{1,2\}$). We define:

$$ \overline{x_i} = rac{x_i - min_i}{max_i - min_i} $$

This normalization maps all cost values into the interval $[0...1]$.

## An Algorithm for Solving BCP-BOS
We now present BCP-BOA*, an algorithm that solves BCP-BOS. BCP-BOA* uses the BO-BFS framework described above but modifies it to fit BCP-BOS.

BCP-BOA* executes a best-first search, which maintains an OPEN list and a CLOSED list. Each node contains a vertex $v$, a g-value $(g_1,g_2)$ and an h-value $(h_1,h_2)$. BCP-BOA* starts by inserting a root node $r$ into OPEN. When a node $n$ that contains vertex $v$ is expanded, BCP-BOA* moves $n$ from OPEN to CLOSED and creates a new node $n^{\prime}$ for each successor vertex $v^{\prime}$ of $v$. 

Next, to address the bounded-cost requirement for each new node $n^{\prime}$, BCP-BOA* discards $n^{\prime}$ if one of its $f_i$-values exceeds the bounds (i.e., if $f_1(n^{\prime}) = g_1(n^{\prime}) + h_1(n^{\prime}) > b_1$ or if $f_2(n^{\prime}) = g_2(n^{\prime}) + h_2(n^{\prime}) > b_2$). This step is called **bound pruning**. If both costs are within their bounds, BCP-BOA* performs a dominance check.

## Choosing the Best Node in OPEN
Given the normalized bounds $\overline{b_1}$ and $\overline{b_2}$, BCP-BOA* has to search the rectangle whose extreme points are $(0,0)$ and $(\overline{b_1},\overline{b_2})$. We describe strategies that search this rectangle systematically using an ordering function $O$.

### Lexicographic Orderings
We define the Lex1 and Lex2 ordering functions as follows:
* $O_{lex1}(n,m) = O(f_1, f_2, n, m)$
* $O_{lex2}(n,m) = O(f_2, f_1, n, m)$

$O_{lex1}$ first expands nodes along the left vertical line from bottom to top, then the second left vertical line from bottom to top, and so on. $O_{lex2}$ first expands nodes along the bottom horizontal line from left to right. When $b_1$ is larger than $b_2$, $O_{lex2}$ often outperforms $O_{lex1}$ and vice versa.

### Minimum and Maximum Orderings
Let $F_{min}(n) = \min(\overline{f}_1(n), \overline{f}_2(n))$ and $F_{max}(n) = \max(\overline{f}_1(n), \overline{f}_2(n))$. We define the Min and Max ordering functions as follows:
* $O_{min}(n,m) = O(F_{min}, F_{max}, n, m)$
* $O_{max}(n,m) = O(F_{max}, F_{min}, n, m)$

### Average Ordering
Let $F_{avg}(n) = (\overline{f}_1(n) + \overline{f}_2(n))/2$. We define the Average ordering function as follows:
* $O_{avg}(n,m) = O(F_{avg}, F_{min}, n, m)$

## Experimental Results
We compared all variants of BCP-BOA* on the BAY road map. We chose five nodes from the POF as pivots to set the two bounds between $min_i$ and $max_i$:
1. **FTL**: farthest top-left POF solution path
2. **FBR**: farthest bottom-right POF solution path
3. **MD**: POF solution path with the smallest difference between its $\overline{c_1}$ and $\overline{c_2}$ costs
4. **TL**: closest to the middle between MD and FTL
5. **BR**: closest to the middle between MD and FBR

**Table 1: Number of expanded nodes (in thousands) for Normalized BCP-BOA***

| Ordering | Z2: FTL | Z2: TL | Z2: MD | Z2: BR | Z2: FBR | Z3: FTL | Z3: TL | Z3: MD | Z3: BR | Z3: FBR | Z4: FTL | Z4: TL | Z4: MD | Z4: BR | Z4: FBR | Z5: Any |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **$b_1$** | .26 | .42 | .51 | .60 | .96 | .50 | .61 | .67 | .73 | .97 | .74 | .80 | .83 | .85 | .98 | 1.00 |
| **$b_2$** | .95 | .60 | .51 | .43 | .21 | .97 | .74 | .76 | .62 | .50 | .98 | .87 | .84 | .81 | .75 | 1.00 |
| **Lex1** | 1.1 | 83.6 | 105.8 | 110.3 | 81.6 | 1.1 | 61.1 | 88.2 | 107.6 | 106.9 | 1.1 | 26.7 | 52.6 | 55.8 | 84.3 | 1.1 |
| **Lex2** | 117.7 | 107.4 | 78.3 | 49.1 | 0.8 | 103.7 | 68.5 | 44.3 | 24.6 | 0.8 | 38.5 | 19.5 | 11.3 | 6.3 | 0.8 | 0.8 |
| **Sel. Lex** | 1.1 | 83.6 | 92.8 | 49.1 | 0.8 | 1.1 | 61.1 | 59.7 | 24.6 | 0.8 | 1.1 | 26.6 | 16.5 | 6.3 | 0.8 | 0.8 |
| **Min** | 1.5 | 96.6 | 107.2 | 65.3 | 0.9 | 1.8 | 60.4 | 71.0 | 36.3 | 1.0 | 2.1 | 18.7 | 15.3 | 9.5 | 1.2 | 1.3 |
| **Max** | 117.8 | 123.8 | 124.5 | 120.5 | 87.9 | 123.7 | 124.0 | 124.5 | 124.2 | 120.6 | 124.5 | 124.5 | 124.5 | 124.5 | 124.5 | 124.4 |
| **Average** | 141.6 | 174.7 | 175.4 | 161.1 | 113.6 | 173.6 | 180.7 | 183.4 | 179.1 | 166.0 | 181.0 | 182.1 | 181.3 | 181.3 | 182.1 | 181.3 |
| **All POF** | 162.8 | 226.0 | 241.5 | 245.9 | 202.8 | 264.4 | 303.9 | 312.4 | 314.0 | 306.1 | 337.6 | 349.2 | 353.6 | 356.3 | 354.0 | 369.4 |
| **#POF Sol.** | 67 | 56 | 56 | 51 | 50 | 98 | 92 | 88 | 84 | 86 | 124 | 120 | 118 | 115 | 116 | 147 |

Selective Lex is an intelligent variant which selects Lex2 if $b_1 > b_2$ and Lex1 otherwise. Thus, it exploits the benefits of both Lex variants. Our results confirm that Selective Lex was the most robust variant: it was either the best ordering function or very close to it.

## Conclusions
We presented BCP-BOA* and several variants of it, and concluded that Selective Lex is the best ordering function for ordering nodes in OPEN. It is future work to lift the requirement that the returned solution path must be in the POF, thereby allowing variants of Potential Search (Stern et al. 2014) to be used. Our normalization scheme can be adapted to other BOS algorithms. Finally, our work can be generalized to Multi-Objective Search.

## Acknowledgements
This research was supported by the United States-Israel Binational Science Foundation (BSF) under grant numbers 2017692 and 2021643, Israel Science Foundation (ISF) under grant number 844/17, the National Science Foundation (NSF) under grant numbers 1409987, 1724392, 1817189, 1837779, 1935712, 2112533, 2121028, and by Centro Nacional de Inteligencia Artificial under grant number FB210017.

## References
* Breugem, T.; Dollevoet, T.; and Heuvel, W. 2017. Analysis of FPTASes for the Multi-Objective Shortest Path Problem. Computers & Operations Research, 78: 44-58.
* Bronfman, A.; Marianov, V.; Paredes-Belmar, G.; and Lüer-Villagra, A. 2015. The Maximin HAZMAT Routing Problem. European Journal of Operational Research, 241(1): 15-27.
* Ehrgott, M. 2005. Multicriteria Optimization, volume 491. Springer Science & Business Media.
* Fu, M.; Kuntz, A.; Salzman, O.; and Alterovitz, R. 2019. Toward Asymptotically-Optimal Inspection Planning Via Efficient Near-Optimal Graph Search. In Robotics: Science and Systems XV.
* Fu, M.; Salzman, O.; and Alterovitz, R. 2021. Computationally-Efficient Roadmap-Based Inspection Planning via Incremental Lazy Search. In ICRA, 7449-7456.
* Goldin, B.; and Salzman, O. 2021. Approximate Bi-Criteria Search by Efficient Representation of Subsets of the Pareto-Optimal Frontier. In ICAPS, 149-158.
* Hart, P. E.; Nilsson, N. J.; and Raphael, B. 1968. A Formal Basis for the Heuristic Determination of Minimum Cost Paths. IEEE Transactions on Systems Science and Cybernetics, 4(2): 100-107.
* Hernández Ulloa, C.; Yeoh, W.; Baier, J. A.; Zhang, H.; Suazo, L.; and Koenig, S. 2020. A Simple and Fast Bi-Objective Search Algorithm. In ICAPS, 143-151.
* Mandow, L.; and De La Cruz, J. L. P. 2005. A New Approach to Multi-Objective A* Search. In IJCAI, 218-223.
* Pulido, F.-J.; Mandow, L.; and Pérez-de-la Cruz, J.-L. 2015. Dimensionality Reduction in Multiobjective Shortest Path Search. Computers & Operations Research, 64: 60-70.
* Raith, A.; and Ehrgott, M. 2009. A Comparison Cf Solution Strategies for Biobjective Shortest Path Problems. Computers & Operations Research, 36(4): 1299-1331.
* Stern, R.; Felner, A.; van den Berg, J.; Puzis, R.; Shah, R.; and Goldberg, K. 2014. Potential-Based Bounded-Cost Search and Anytime Non-Parametric A*. Artificial Intelligence, 214: 1-25.
* Stewart, B. S.; and White III, C. C. 1991. Multiobjective A*. Journal of the ACM, 38(4): 775-814.
