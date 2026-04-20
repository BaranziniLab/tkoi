# Gene Metadata

A dataset containing information about genes, including identifiers,
names, and degree in the network.

## Usage

``` r
genes
```

## Format

A data frame with the following columns:

- id:

  Character. Unique identifier for the gene in the network.

- identifier:

  Integer. Internal identifier for the gene.

- ensembl:

  Character. Ensembl gene identifier.

- name:

  Character. Gene symbol or common name.

- degree:

  Integer. Degree of the gene node in the network, representing the
  number of connections.

## Details

This dataset is essential for integrating gene-level information with
biological networks and functional studies.

## Examples

``` r
data(genes)
head(genes)
#>                                                id identifier         ensembl
#>                                            <char>      <int>          <char>
#> 1:   4:c77f6410-bc08-43ba-a172-0503ab1c93db:16050     389611 ENSG00000233132
#> 2:   4:c77f6410-bc08-43ba-a172-0503ab1c93db:17633     389630 ENSG00000230045
#> 3:   4:c77f6410-bc08-43ba-a172-0503ab1c93db:47899     151507 ENSG00000293137
#> 4: 4:c77f6410-bc08-43ba-a172-0503ab1c93db:2490276          1 ENSG00000121410
#> 5: 4:c77f6410-bc08-43ba-a172-0503ab1c93db:2490277          2 ENSG00000175899
#> 6: 4:c77f6410-bc08-43ba-a172-0503ab1c93db:2490278          9 ENSG00000171428
#>        name degree
#>      <char>  <int>
#> 1:  FAM90A3      8
#> 2: FAM90A15     26
#> 3:    MSL3B    105
#> 4:     A1BG    129
#> 5:      A2M    566
#> 6:     NAT1    239
```
