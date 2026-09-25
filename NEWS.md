# tkoi 1.1.0

## Performance

* The core of `run_tkoi()` has been rewritten in C++ (via Rcpp). Personalized
  PageRank for the observed data and for every permutation is solved in
  batches of up to 16 vectors per pass over the network, with conjugate
  gradient on multiple threads. Two-hop neighbourhood counts use bitsets, and
  the degree-matched null gene sets are sampled in C++.
* On an Apple M2 (8 cores), with a 21,414-gene differential expression
  dataset (1,201 seed genes) on the full `tkoi_net`:
  * 10 permutations take 5.3 s (48.5 s with tkoi 1.0.0). The first run in a
    session, which also builds the network matrix, takes 8.8 s, and a run on
    one core (`n_cores = 1`) takes 10.9 s.
  * 100 permutations take 25.5 to 38.1 s, depending on thermal throttling
    (244.3 s with tkoi 1.0.0).
  * Overall, runs are about 6 to 10 times faster, and a single core is still
    faster than tkoi 1.0.0.
* The statistics match tkoi 1.0.0 run with the same seed: the largest
  difference in `beta` is 7.6e-9 with 10 permutations and 5.4e-9 with 100
  permutations, and PageRank vectors agree to within 6e-15. The result tables
  differ from 1.0.0 only through the deliberate changes listed under
  "Changes to result tables"; apart from those, rows can be ordered
  differently only where FDR and `beta` agree to 10 significant digits.
* The network matrix is built once and cached, so repeated runs on the same
  network in a session skip this step. The cache is keyed on the graph
  structure (`igraph::graph_id()`) and the `weight` edge attribute; vertex
  metadata (names, types, identifiers, degrees) is read afresh on every call,
  so edits to vertex attributes are always honoured.

## Changes to result tables

These changes make the reported statistics more complete and more accurate,
and change `network_summary_statistics` compared with tkoi 1.0.0:

* Every node is reported. tkoi 1.0.0 inner-joined each node-type table with
  its annotation dataset, silently dropping nodes without a curated
  annotation (about 20% of Disease, 71% of MiRNA, and 10% of Gene nodes on
  `tkoi_net`), even though they were scored. Tables are now left-joined, so
  these nodes appear with missing annotation columns (and one row per node:
  duplicated annotation identifiers no longer duplicate rows).
* The Compound FDR is adjusted over the reported compounds only. tkoi 1.0.0
  adjusted p-values over all 554,526 Compound nodes but reported only human
  metabolites (about 24,000), making the reported FDRs about 20 times too
  conservative. Compounds that are not human metabolites are now left out
  before the adjustment.
* Nodes whose null PageRank is constant (standard deviation 0) are
  untestable: their `beta`, `p_value`, and `fdr` are `NaN`, they are left out
  of the FDR adjustment, and they are listed last. In tkoi 1.0.0 a seed gene in
  a small disconnected component that no replacement gene reaches got
  `beta = Inf`, `p_value = 0`, and ranked first. `compute_network_enrichment()`
  follows the same rule.
* `human_metabolites` was repaired. Its ChEBI entries had a doubled prefix
  (`"CHEBI:CHEBI:50599"`), its ChEMBL entries lacked the `CHEMBL` prefix, and
  it contained placeholders (`"CHEBI:NA"`, a bare `"inchikey:"`) and
  duplicates, so only InChIKey-identified compounds could pass the human
  metabolite filter. The vector now holds 231,466 unique identifiers in the
  form of Compound node identifiers, and 45 ChEBI-identified compounds are
  reported in addition to the 24,099 InChIKey-identified ones.
* `visualize_topn()` labels nodes without a name by their identifier.

## New arguments

* `run_tkoi()` gains:
  * `n_cores = NULL`: number of threads. `NULL` uses every available core, or
    `getOption("tkoi.n_cores")` when that option is set. Requests above the
    number of available cores are capped; the count of available cores
    respects Slurm (`SLURM_CPUS_PER_TASK`) and Linux cgroup CPU quotas.
    Threads are used on every platform, including Windows. (Development
    builds of 1.0.0 had an `n_cores` argument that defaulted to `1` and used
    `parallel::mclapply()`, which worked only on macOS and Linux and gave
    results that depended on the number of cores.)
  * `keep_permutations = TRUE`: set to `FALSE` to store only the null mean and
    standard deviation (`null_mean`, `null_sd`) in `pagerank_data` instead of
    every permutation, which saves memory for large `n_permutation`.
  * `tolerance = 1e-14`: relative residual at which the PageRank solver
    stops. A vector is accepted when both the 2-norm residual of the symmetric
    system the solver works on and the 1-norm residual of the original
    PageRank system (relative to the personalization vector) are at most
    `tolerance`; the second condition bounds the 1-norm error of the PageRank
    vector by about `2 * tolerance`. The solver needs about 55 iterations per
    vector on `tkoi_net`. At this default, against a reference solution
    computed at `1e-16` on `tkoi_net`, the PageRank vectors are at least as
    accurate as those of `igraph::page_rank()`, per node as well as overall
    (1-norm error 2.3e-15 versus 4.1e-12, and largest per-node relative error
    7.5e-8 versus 2.7e-7).
  * `verbose = TRUE`: set to `FALSE` to silence progress messages.
* `plot_network()` gains a `subnetwork` argument (default `tkoi::tkoi_net`),
  so it can be used with the network that was passed to `run_tkoi()`.

## Breaking and behaviour changes

* `tkoi_net` is now shipped as a plain igraph object. The network is no longer
  encrypted: the `sodium` dependency, `R/sysdata.rda`, and the
  `network_attributes` dataset have been removed, and passing a raw
  (encrypted) vector as `subnetwork` is an error. Loading `tkoi_net` takes a
  few seconds and about 0.5 GB of memory.
* `run_tkoi()` now uses every available core by default. Set `n_cores = 1`
  or `options(tkoi.n_cores = 1)` to use a single thread.
* Reproducibility: `set.seed()` fully determines the result. Results are
  bit-identical for any `n_cores`, and do not depend on how many PageRank
  vectors are solved together (a block width chosen from the available
  memory). On `tkoi_net`, the same seed draws the same null seed sets as a
  sequential tkoi 1.0.0 run.
* Result tables order nodes whose FDR and `beta` agree to 10 significant
  digits deterministically (in network order), so row order no longer depends
  on floating-point noise.
* `maximum_iteration` is now honoured. It was previously ignored by igraph's
  PRPACK solver. A warning is raised if a PageRank vector does not reach
  `tolerance` within `maximum_iteration` iterations.
* `maximum_iteration` and `n_permutation` must not exceed
  `.Machine$integer.max`.
* `run_tkoi()` checks the memory a run needs before the null seed sets are
  drawn and before any PageRank is computed. The check uses physical memory,
  or a Linux cgroup memory limit when that is lower, and stops with advice
  (for example `keep_permutations = FALSE`) if the run would not fit.
* `pagerank_data` now has node IDs as row names.
* For a custom `subnetwork`, replacement genes in the permutation null are
  drawn only from genes that are vertices of that network. For `tkoi_net`
  this is every gene in `genes`, so results are unchanged.
* Nodes with a missing or empty type are grouped under the type `Unknown`.
  Node types without a curated annotation (possible with a custom
  `subnetwork`) are kept in `network_summary_statistics`, without annotation
  columns.
* R >= 4.1 is required.

## Bug fixes

* For a custom `subnetwork`, seed genes that are not vertices of the network
  are now dropped (tkoi 1.0.0 stopped with an error). They are dropped before
  the null is sampled, so the observed run and every null run use the same
  seeds with the same weights and the null is not biased.
* Rows with a missing or blank `gene_name`, and repeated rows for the same
  gene, are now handled the same way by `run_tkoi()`,
  `export_gene_exploration_data()`, `make_gene_exploration_plot()`, and
  `plot_network()`: blank rows are ignored (they previously matched every
  gene without an Ensembl ID) and only the first row of each gene is used.
* `plot_network()` works again. It previously failed because it called
  `get_neighboring_nodes()` without a network and read the encrypted
  `tkoi_net` directly. Paths of up to two hops are now found from neighbour
  sets, the log fold change color scale is centred on zero, and a clear error
  is given when no significant gene is within reach of the target. On custom
  networks it uses vertex names when there is no `identifier` attribute and
  reads node types from `labels` or `label`.
* `get_neighboring_nodes()` now defaults to `subnetwork = tkoi::tkoi_net`
  (which works now that the network is not encrypted) and validates its
  inputs.
* `run_gene_enrichment()` no longer runs the unused
  `enrichplot::pairwise_termsim()` step and handles GO ontologies that return
  no enrichment results.
* `make_gene_exploration_plot()` no longer prints a ggplot2 message about an
  unknown label.
* `visualize_topn()` validates `category` and `top_n` and tolerates duplicated
  labels.
* `compute_network_enrichment()` also accepts a data frame with one row per
  node, such as the `pagerank_data` slot of a run with
  `keep_permutations = TRUE`, or a list of equal-length vectors, which is
  treated the same way. It gives a clear error when `node` has no `pagerank`
  value or fewer than two `perm*` values.
* `run_tkoi()` gives clear errors for invalid inputs, such as missing or
  non-numeric columns, out-of-range thresholds, no genes passing the filters,
  or infinite `logfc` values among the selected genes.

## Documentation

* The `genes`, `human_metabolites`, `compound_annotation`, and `tkoi_net`
  help pages describe the identifiers the datasets actually contain, and the
  `export_gene_exploration_data()`, `make_gene_exploration_plot()`, and
  `export_network_summary_statistics()` help pages have examples based on a
  `run_tkoi()` result instead of hand-built objects.

## Internal

* New C++ engine in `src/tkoi_engine.cpp`, with R-side helpers in
  `R/engine.R`. The package now needs a compiler toolchain to install from
  source (Xcode Command Line Tools on macOS, Rtools on Windows).
* The package lives at <https://github.com/BaranziniLab/tkoi>; the former
  `Broccolito/tkoi` address redirects there. Bioconductor dependencies are
  resolved automatically through `biocViews`.
* Dependencies `sodium`, `purrr`, `readxl`, and `enrichplot` were removed;
  `data.table` moved to Suggests; `Rcpp` was added.
* Threads are capped at 2 when `_R_CHECK_LIMIT_CORES_` is set, as on CRAN.
* New test suite built on a small fixture network, including parity tests
  against the tkoi 1.0.0 algorithm. Tests that need the full `tkoi_net` or
  `clusterProfiler` run only when `TKOI_EXTENDED_TESTS=true`.

# tkoi 1.0.0

* Initial release.
