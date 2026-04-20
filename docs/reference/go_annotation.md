# Gene Ontology (GO) Term Annotations

A comprehensive annotation table of Gene Ontology (GO) terms, including
biological process, molecular function, and cellular component
namespaces. This dataset supports the interpretation and functional
enrichment of graph nodes labeled with GO identifiers within the tKOI
framework.

## Usage

``` r
go_annotation
```

## Format

A data frame with 47,921 rows and 4 columns:

- identifier:

  A character string representing the GO identifier (e.g.,
  "GO:0000001").

- name:

  The name of the GO term (e.g., "mitochondrion inheritance").

- namespace:

  The ontology category to which the term belongs, one of
  `"biological_process"`, `"molecular_function"`, or
  `"cellular_component"`.

- definition:

  A textual definition of the GO term, often including references and
  source attributions.

## Source

Gene Ontology Consortium <http://geneontology.org/>

## Details

The dataset enables mapping of functional annotations to nodes in a
biological knowledge graph and is particularly useful for characterizing
enriched terms from network-based analyses such as those performed with
[`run_tkoi`](run_tkoi.md). GO terms marked as "obsolete" are included
for completeness but may require special handling during analysis.

## See also

[`anatomy_annotation`](anatomy_annotation.md),
[`celltype_annotation`](celltype_annotation.md),
[`run_tkoi`](run_tkoi.md)

## Examples

``` r
data(go_annotation)
table(go_annotation$namespace)
#> 
#> biological_process cellular_component           external molecular_function 
#>              30741               4525                  1              12654 
head(go_annotation[go_annotation$namespace == "molecular_function", ])
#>    identifier                                                     name
#>        <char>                                                   <char>
#> 1: GO:0000005                    obsolete ribosomal chaperone activity
#> 2: GO:0000006    high-affinity zinc transmembrane transporter activity
#> 3: GO:0000007 low-affinity zinc ion transmembrane transporter activity
#> 4: GO:0000008                                     obsolete thioredoxin
#> 5: GO:0000009                   alpha-1,6-mannosyltransferase activity
#> 6: GO:0000010                heptaprenyl diphosphate synthase activity
#>             namespace
#>                <char>
#> 1: molecular_function
#> 2: molecular_function
#> 3: molecular_function
#> 4: molecular_function
#> 5: molecular_function
#> 6: molecular_function
#>                                                                                                                                                                                                                                                                                                                                                                                                       definition
#>                                                                                                                                                                                                                                                                                                                                                                                                           <char>
#> 1:                                                                                                                                                                                              OBSOLETE. Assists in the correct assembly of ribosomes or ribosomal subunits in vivo, but is not a component of the assembled ribosome when performing its normal biological function."" [GOC:jl, PMID:12150913]
#> 2:                                                                                                                                     Enables the transfer of zinc ions (Zn2+) from one side of a membrane to the other, probably powered by proton motive force. In high-affinity transport the transporter is able to bind the solute even if it is only present at very low concentrations."" [TC:2.A.5.1.1]
#> 3:                                                                          Enables the transfer of a solute or solutes from one side of a membrane to the other according to the reaction: Zn2+ = Zn2+, probably powered by proton motive force. In low-affinity transport the transporter is able to bind the solute only if it is present at very high concentrations."" [GOC:mtg_transport, ISBN:0815340729]
#> 4: OBSOLETE. A small disulfide-containing redox protein that serves as a general protein disulfide oxidoreductase. Interacts with a broad range of proteins by a redox mechanism, based on the reversible oxidation of 2 cysteine thiol groups to a disulfide, accompanied by the transfer of 2 electrons and 2 protons. The net result is the covalent interconversion of a disulfide and a dithiol."" [GOC:kd]
#> 5:                                                                                                                                                                                                                                                                              Catalysis of the transfer of a mannose residue to an oligosaccharide, forming an alpha-(1->6) linkage."" [GOC:mcc, PMID:2644248]
#> 6:                                                                                                                                                                                                                                         Catalysis of the reaction: (2E,6E)-farnesyl diphosphate + 4 isopentenyl diphosphate = 4 diphosphate + all-trans-heptaprenyl diphosphate."" [PMID:9708911, RHEA:27794]
```
