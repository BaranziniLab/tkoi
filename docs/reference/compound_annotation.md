# Chemical Compound Annotations

A curated annotation table for small molecules and chemical compounds
used in pharmacological, biochemical, or research contexts. Each entry
includes a standardized identifier (e.g., InChIKey or ChEMBL ID) and the
corresponding compound name, when available. This dataset is used within
the tKOI framework to annotate nodes representing chemical entities or
drugs in biological knowledge graphs.

## Usage

``` r
compound_annotation
```

## Format

A data frame with 554526 rows and 2 columns:

- identifier:

  A character string representing the unique compound ID, either as an
  InChIKey (e.g., "inchikey:NZLJDTKLZIMONR-UHFFFAOYSA-N") or a ChEMBL
  compound accession (e.g., "chembl.compound:CHEMBL5219790").

- name:

  The human-readable compound name, when available (e.g.,
  "ACRIFLAVINE"). Missing values indicate uncharacterized or unnamed
  entries.

## Source

Data compiled from chemical databases such as ChEMBL and InChI registry.

## Details

This dataset enables mapping of compound-level features in biomedical
networks, including chemical perturbagens, drug candidates, or
environmental exposures. It can be joined to tKOI network results using
the `identifier` field to enrich nodes of type "Compound" with
interpretable names.

## See also

[`complex_annotation`](complex_annotation.md),
[`clinicallab_annotation`](clinicallab_annotation.md),
[`run_tkoi`](run_tkoi.md)

## Examples

``` r
data(compound_annotation)
subset(compound_annotation, grepl("CETRIMIDE", name))
#>                      identifier      name
#> 3 chembl.compound:CHEMBL5219790 CETRIMIDE
#> 4 chembl.compound:CHEMBL5219942 CETRIMIDE
```
