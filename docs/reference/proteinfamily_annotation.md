# Protein Family Annotations (Pfam Clans)

A curated annotation table of protein families, primarily derived from
Pfam clan classifications. Each entry maps a protein family identifier
to a short descriptive name, supporting the annotation of graph nodes
that represent protein families in the tKOI knowledge network.

## Usage

``` r
proteinfamily_annotation
```

## Format

A data frame with 659 rows and 2 columns:

- identifier:

  A character string representing the Pfam clan or protein family
  identifier (e.g., "CL0163").

- name:

  The short name of the protein family (e.g., "Calcineurin"). These
  often reflect structural, evolutionary, or functional groupings of
  related protein domains.

## Source

Pfam Clans <https://pfam.xfam.org/clan>

## Details

Protein families group related domains based on common evolutionary
origins, structure, or function. This dataset allows high-level
functional grouping of nodes in biological networks and supports
annotation of gene products in systems biology and bioinformatics
analyses.

The `identifier` field serves as a join key for nodes of type
"ProteinFamily" in the tKOI framework.

## See also

[`proteindomain_annotation`](proteindomain_annotation.md),
[`protein_annotation`](protein_annotation.md), [`run_tkoi`](run_tkoi.md)

## Examples

``` r
data(proteinfamily_annotation)
subset(proteinfamily_annotation, grepl("Calcineurin", name))
#>   identifier        name
#> 2     CL0163 Calcineurin
```
