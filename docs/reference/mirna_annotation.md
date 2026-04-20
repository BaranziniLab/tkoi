# Human miRNA Annotations (miRBase-style)

A curated dataset mapping human mature miRNA identifiers to their
corresponding miRBase precursor accessions.

## Usage

``` r
mirna_annotation
```

## Format

A data frame with two columns and multiple rows (1 per miRNA):

- identifier:

  Character string representing the mature miRNA ID, including the
  species prefix (e.g., `"hsa-miR-6129"`).

- accession:

  miRBase accession number for the precursor miRNA hairpin (e.g.,
  `"MI0021274"`).

## Source

Extracted from a parsed version of miRBase release data.

## Details

This dataset enables annotation and cross-referencing of mature miRNA
names with miRBase precursor IDs. It is useful for linking results from
miRNA expression studies with external reference databases such as
miRBase.

The `identifier` field typically corresponds to results from expression
matrices, differential analysis tools, or network analysis pipelines.

## See also

[`run_tkoi`](run_tkoi.md),
[`reaction_annotation`](reaction_annotation.md)

## Examples

``` r
data(mirna_annotation)
head(mirna_annotation)
#>     identifier accession
#> 1 hsa-miR-6129 MI0021274
#> 2 hsa-miR-6133 MI0021278
#> 3 hsa-miR-4510 MI0016876
#> 4  hsa-miR-23c MI0016010
#> 5 hsa-miR-4292 MI0015897
#> 6 hsa-miR-6127 MI0021271
subset(mirna_annotation, grepl("miR-23", identifier))
#>       identifier accession
#> 4    hsa-miR-23c MI0016010
#> 479 hsa-miR-2392 MI0016870
```
