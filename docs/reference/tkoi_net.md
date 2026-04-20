# Human-Specific Heterogeneous Network

An igraph object representing a human-specific heterogeneous biological
network.

## Usage

``` r
tkoi_net
```

## Format

An igraph object where nodes and edges represent biological entities and
their interactions. The node types included are:

- Anatomy:

  Nodes representing anatomical structures and systems.

- BiologicalProcess:

  Nodes for functional biological processes, such as signaling pathways.

- CellType:

  Nodes describing different cell types.

- CellularComponent:

  Nodes for subcellular structures, organelles, and macromolecular
  complexes.

- ClinicalLab:

  Nodes representing clinical measurements and diagnostic data.

- Complex:

  Nodes for molecular and protein complexes.

- Compound:

  Nodes for endogenous metabolites in human.

- Disease:

  Nodes for diseases and pathological conditions.

- EC:

  Nodes categorized by Enzyme Commission numbers.

- Gene:

  Nodes for genetic elements, such as genes and genetic markers.

- MiRNA:

  Nodes for microRNAs and their regulatory roles.

- MolecularFunction:

  Nodes describing molecular activities performed by gene products.

- Pathway:

  Nodes representing sequences of molecular interactions and reactions.

- Protein:

  Nodes for protein molecules.

- ProteinDomain:

  Nodes for specific structural or functional domains within proteins.

- ProteinFamily:

  Nodes for groups of evolutionarily related proteins.

- PwGroup:

  Nodes for pathway groups aggregating multiple related pathways.

- Reaction:

  Nodes for biochemical reactions and their participants.

## Details

This heterogeneous network integrates multiple biological datasets to
represent complex relationships within the human system. It serves as
the foundation for network-based analyses in the `tkoi` package, such as
personalized PageRank calculations and enrichment analyses.

## Examples

``` r
data(tkoi_net)
igraph::summary(tkoi_net)
#> Error: 'summary' is not an exported object from 'namespace:igraph'
igraph::V(tkoi_net)$name[1:10]
#> Error in ensure_igraph(graph): Must provide a graph object (provided wrong object type).
igraph::E(tkoi_net)[1:10]
#> Error in ensure_igraph(graph): Must provide a graph object (provided wrong object type).
```
