# Anatomical Entity Annotations from UBERON

A curated annotation table of anatomical entities sourced from the
UBERON ontology. This dataset provides standard identifiers,
human-readable names, and textual definitions for anatomical structures.
It is used within the tKOI analysis pipeline to enrich network nodes of
type "Anatomy" with biological context and ontology-based metadata.

## Usage

``` r
anatomy_annotation
```

## Format

A data frame with 15,668 rows and 3 columns:

- identifier:

  A character string representing the UBERON ontology identifier (e.g.,
  "UBERON:0000002").

- name:

  The human-readable name of the anatomical structure (e.g., "uterine
  cervix").

- definition:

  A character string describing the biological definition or functional
  role of the anatomical entity.

## Source

UBERON Ontology <http://uberon.github.io/>

## Details

The dataset supports the annotation and interpretation of graph nodes
labeled as anatomical concepts in biological networks. It can be joined
to results using the `identifier` column to enhance interpretability of
enriched or prioritized nodes.

## See also

[`run_tkoi`](run_tkoi.md),
[`celltype_annotation`](celltype_annotation.md),
[`disease_annotation`](disease_annotation.md)

## Examples

``` r
data(anatomy_annotation)
head(anatomy_annotation)
#>        identifier                name
#>            <char>              <char>
#> 1: UBERON:0000000   processual entity
#> 2: UBERON:0000002      uterine cervix
#> 3: UBERON:0000003               naris
#> 4: UBERON:0000004                nose
#> 5: UBERON:0000005  chemosensory organ
#> 6: UBERON:0000006 islet of Langerhans
#>                                                                                                                                         definition
#>                                                                                                                                             <char>
#> 1: An occurrent [span:Occurrent] that exists in time by occurring or happening, has temporal parts and always involves and depends on some entity.
#> 2:                                                              Lower, narrow portion of the uterus where it joins with the top end of the vagina.
#> 3:                                      Orifice of the olfactory system. The naris is the route by which odorants enter the olfactory system[MAH].
#> 4:     The olfactory organ of vertebrates, consisting of nares, olfactory epithelia and the structures and skeletal framework of the nasal cavity.
#> 5:                                                                                                                                                
#> 6:                                                             The clusters of hormone-producing cells that are scattered throughout the pancreas.
```
