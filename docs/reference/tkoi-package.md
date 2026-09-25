# tkoi: Transcriptomic Knowledge-Graph Omics Integration

Network-aware gene enrichment analysis. Transcriptomic signals are
propagated through a human biological knowledge graph with personalized
PageRank, and a degree-matched permutation null identifies enriched
biological concepts. The main entry point is
[`run_tkoi()`](run_tkoi.md).

## Performance

The PageRank solver, neighborhood counting, and null sampling are
implemented in C++ and run on multiple threads. See the `n_cores`
argument of [`run_tkoi()`](run_tkoi.md) and the `tkoi.n_cores` option.

## See also

Useful links:

- <https://github.com/BaranziniLab/tkoi>

- Report bugs at <https://github.com/BaranziniLab/tkoi/issues>

## Author

**Maintainer**: Wanjun Gu <wanjun.gu@ucsf.edu>
([ORCID](https://orcid.org/0000-0002-7342-7000))
