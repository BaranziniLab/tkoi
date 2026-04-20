# Disease Ontology (DO) Annotations

A structured annotation table mapping Disease Ontology (DOID)
identifiers to human-readable disease names and definitions. This
dataset is used in the tKOI framework to annotate network nodes
associated with disease terms, enabling interpretation of enriched or
prioritized disease-related features in biological networks.

## Usage

``` r
disease_annotation
```

## Format

A data frame with 10002 rows and 3 columns:

- identifier:

  A character string representing the Disease Ontology ID (e.g.,
  "DOID:0040003").

- name:

  The name of the disease (e.g., "benzylpenicillin allergy").

- definition:

  A detailed textual definition of the disease, including etiology,
  affected systems, or allergy triggers, often with references to
  biomedical sources.

## Source

Disease Ontology <https://disease-ontology.org/>

## Details

Disease Ontology (DO) provides a standardized, ontology-based
representation of human diseases. This dataset allows integration of
disease concepts into biological graphs and supports downstream analysis
of disease associations.

The `identifier` column serves as a join key to nodes of type "Disease"
in tKOI networks, while the `definition` column provides biological and
clinical context.

## See also

[`clinicallab_annotation`](clinicallab_annotation.md),
[`compound_annotation`](compound_annotation.md),
[`run_tkoi`](run_tkoi.md)

## Examples

``` r
data(disease_annotation)
subset(disease_annotation, grepl("allergy", name))
#>        identifier                     name
#>            <char>                   <char>
#>   1: DOID:0040001           shrimp allergy
#>   2: DOID:0040002          aspirin allergy
#>   3: DOID:0040003 benzylpenicillin allergy
#>   4: DOID:0040004      amoxicillin allergy
#>   5: DOID:0040005      ceftriaxone allergy
#>  ---                                      
#>  99:    DOID:3660            wheat allergy
#> 100:    DOID:4376             milk allergy
#> 101:    DOID:4377              egg allergy
#> 102:    DOID:4378           peanut allergy
#> 103:    DOID:4379              nut allergy
#>                                                                                                                                                                                                                                                                                                                                                                      definition
#>                                                                                                                                                                                                                                                                                                                                                                          <char>
#>   1:                                                                                                                                                                                                                                      ""A crustacean allergy that has_allergic_trigger shrimp."" [url:https\\://www.ncbi.nlm.nih.gov/pubmed/20471069] {comment=""IEDB:RV""}
#>   2:                                                                                                                                                                                                                               ""A drug allergy that has_allergic_trigger acetylsalicylic acid."" [url:https\\://www.ncbi.nlm.nih.gov/pubmed/2468301] {comment=""IEDB:RV""}
#>   3:                                                                                                                                                                                                                           ""A beta-lactam allergy that has_allergic_trigger benzylpenicillin."" [url:https\\://www.ncbi.nlm.nih.gov/pubmed/14483916] {comment=""IEDB:RV""}
#>   4:                                                                                                                                                                                                                                ""A beta-lactam allergy that has_allergic_trigger amoxicillin."" [url:https\\://www.ncbi.nlm.nih.gov/pubmed/11746950] {comment=""IEDB:RV""}
#>   5:                                                                                                                                                                                                                              ""A cephalosporin allergy that has_allergic_trigger ceftriaxone."" [url:https\\://www.ncbi.nlm.nih.gov/pubmed/12833570] {comment=""IEDB:RV""}
#>  ---                                                                                                                                                                                                                                                                                                                                                                           
#>  99: ""A food allergy that develops_from exposure to and particularly consumption of wheat, and has_symptom that are both gastrointestinal and nongastrointestinal in nature, such as diarrhea, mouth and throat irritation, headache, hives, skin rashes, and anaphylaxis."" [url:https\\://www.mayoclinic.org/diseases-conditions/wheat-allergy/symptoms-causes/syc-20378897]
#> 100:                                                                                                           ""A food allergy that results in adverse immune reaction to one or more of the proteins in cow's milk and/or the milk of other animals, which are normally harmless to the non-allergic individual."" [url:http\\://en.wikipedia.org/wiki/Milk_hypersensitivity]
#> 101:                                                                                         ""A food allergy that is an allergy or hypersensitivity to dietary substances from the yolk or whites of eggs, causing an overreaction of the immune system which may lead to severe physical symptoms."" [url:http\\://en.wikipedia.org/wiki/Allergy#Foods] {comment=""ls:IEDB""}
#> 102:                                                                           ""A legume allergy that is an allergy or hypersensitivity to dietary substances from peanuts causing an overreaction of the immune system which in a small percentage of people may lead to severe physical symptoms."" [url:http\\://en.wikipedia.org/wiki/Allergy#Foods] {comment=""ls:IEDB""}
#> 103:                                                                                                                                            ""A food allergy that develops_from exposure to and particularly consumption of nuts, and has_symptom asthma, skin rashes, throat and eye irritation, and anaphylaxis."" [url:https\\://en.wikipedia.org/wiki/Tree_nut_allergy]
```
