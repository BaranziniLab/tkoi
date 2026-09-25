# Human-Specific Heterogeneous Network

The tKOI knowledge graph: an undirected igraph object with 939,059 nodes
and 10,622,200 edges linking human genes to the concepts they relate to.
It is shipped as a plain igraph object and loads lazily the first time
`tkoi::tkoi_net` is used (this takes a few seconds and about 0.5 GB of
memory).

## Usage

``` r
tkoi_net
```

## Format

An igraph object. Vertex attributes:

- name:

  Unique node ID, e.g. `"4:c77f6410-...:16050"`.

- identifier:

  Source identifier (Entrez ID, GO ID, UBERON ID, ...).

- source:

  Source database.

- labels:

  Node type in Neo4j label form, e.g. `"['Gene']"`.

- degree:

  Node degree in the full network.

The edge attribute `edge_type` names the relationship (e.g.
`"PARTICIPATES_GpBP"`). The node types included are:

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

  Nodes for chemical compounds, identified by InChIKey, ChEBI ID, or
  ChEMBL ID (see [`compound_annotation`](compound_annotation.md)).
  [`run_tkoi`](run_tkoi.md) reports only those in
  [`human_metabolites`](human_metabolites.md).

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
# \donttest{
igraph::vcount(tkoi::tkoi_net)
#> [1] 939059
table(igraph::V(tkoi::tkoi_net)$labels)
#> 
#>           ['Anatomy'] ['BiologicalProcess']          ['CellType'] 
#>                 13770                 12996                  2744 
#> ['CellularComponent']       ['ClinicalLab']           ['Complex'] 
#>                  1708                 59296                  3318 
#>          ['Compound']           ['Disease']                ['EC'] 
#>                554526                 11448                  8764 
#>              ['Gene']             ['MiRNA'] ['MolecularFunction'] 
#>                 19503                  2656                  3569 
#>           ['Pathway']           ['Protein']     ['ProteinDomain'] 
#>                  4831                194076                 14193 
#>     ['ProteinFamily']           ['PwGroup']          ['Reaction'] 
#>                   659                  6343                 24659 
head(igraph::V(tkoi::tkoi_net)$identifier)
#> [1] "UBERON:0003233" "UBERON:2001901" "UBERON:0004321" "UBERON:0002414"
#> [5] "UBERON:2005118" "UBERON:0034769"
# }
```
