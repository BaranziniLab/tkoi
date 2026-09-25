# Chemical Compound Annotations

A curated annotation table for small molecules and chemical compounds
used in pharmacological, biochemical, or research contexts. It has one
row per Compound node of [`tkoi_net`](tkoi_net.md), with the node's
identifier (an InChIKey, a ChEBI ID, or a ChEMBL ID) and the compound
name, when available. This dataset is used within the tKOI framework to
annotate nodes representing chemical entities or drugs in biological
knowledge graphs.

## Usage

``` r
compound_annotation
```

## Format

A data frame with 554,526 rows and 2 columns:

- identifier:

  Character. The compound ID, in one of three forms: an InChIKey
  (530,612 rows, e.g. `"inchikey:NZLJDTKLZIMONR-UHFFFAOYSA-N"`), a ChEBI
  ID (22,652 rows, e.g. `"CHEBI:85476"`), or a ChEMBL compound accession
  (1,262 rows, e.g. `"chembl.compound:CHEMBL5219790"`).

- name:

  The human-readable compound name, when available (e.g.,
  "ACRIFLAVINE"). Missing values indicate uncharacterized or unnamed
  entries.

## Source

Data compiled from chemical databases such as ChEMBL, ChEBI, and the
InChI registry.

## Details

This dataset enables mapping of compound-level features in biomedical
networks, including chemical perturbagens, drug candidates, or
environmental exposures. It can be joined to tKOI network results using
the `identifier` field to enrich nodes of type "Compound" with
interpretable names.

[`run_tkoi`](run_tkoi.md) reports only Compound nodes whose identifier
is in [`human_metabolites`](human_metabolites.md). That vector matches
Compound nodes only through InChIKeys, so the `Compound` table of a
[`run_tkoi()`](run_tkoi.md) result lists InChIKey compounds only; the
ChEBI and ChEMBL rows of this table annotate Compound nodes that
[`run_tkoi()`](run_tkoi.md) does not report.

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

# Identifier forms
table(sub(":.*", "", compound_annotation$identifier))
#> 
#>           CHEBI chembl.compound        inchikey 
#>           22652            1262          530612 
head(subset(compound_annotation, startsWith(identifier, "CHEBI:")))
#>         identifier                            name
#> 527390 CHEBI:85476   O-hydroxyvaleroyl-L-carnitine
#> 527391 CHEBI:48679                   sorbopyranose
#> 527392 CHEBI:60592 poly(fluorene-2,7-diyl) polymer
#> 527393 CHEBI:46687                     diazaalkane
#> 527394 CHEBI:53104                            <NA>
#> 527395 CHEBI:64840                            <NA>
```
