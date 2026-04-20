# Cell Type Ontology Annotations (CL)

A structured annotation table of cell types based on the Cell Ontology
(CL). This dataset provides standardized identifiers, names, and
detailed biological definitions for various cell types, including
epithelial cells, immune cells, and specialized sensory cells. It is
used in the tKOI framework to annotate nodes labeled as cell types and
to support network enrichment and interpretation of cell-specific
biological functions.

## Usage

``` r
celltype_annotation
```

## Format

A data frame with 3,058 rows and 3 columns:

- identifier:

  A character string representing the CL (Cell Ontology) identifier
  (e.g., "CL:0000001").

- name:

  The standardized name of the cell type (e.g., "primary cultured
  cell").

- definition:

  A textual definition of the cell type, including functional roles,
  structural features, and biological context. Obsolete terms are
  clearly indicated.

## Source

Cell Ontology (CL) <https://obophenotype.github.io/cell-ontology/>

## Details

The dataset includes a wide variety of cell types, from generic terms
like "cell" to highly specialized cells such as "tuft cell of olfactory
epithelium." Obsolete terms are retained for ontology completeness but
should be excluded in downstream functional analyses unless specifically
required.

The annotations can be joined to graph-based analysis results using the
`identifier` column, facilitating biological interpretation in network
analysis.

## See also

[`anatomy_annotation`](anatomy_annotation.md),
[`go_annotation`](go_annotation.md), [`run_tkoi`](run_tkoi.md)

## Examples

``` r
data(celltype_annotation)
subset(celltype_annotation, grepl("olfactory", name))
#>     identifier
#>         <char>
#>  1: CL:0000207
#>  2: CL:0000626
#>  3: CL:0000847
#>  4: CL:0000848
#>  5: CL:0000849
#>  6: CL:0000853
#>  7: CL:0002167
#>  8: CL:0002169
#>  9: CL:0002171
#> 10: CL:0002184
#> 11: CL:0011028
#> 12: CL:1001434
#> 13: CL:1001503
#> 14: CL:4030051
#> 15: CL:4030052
#> 16: CL:4052037
#>                                                                      name
#>                                                                    <char>
#>  1:                                               olfactory receptor cell
#>  2:                                                olfactory granule cell
#>  3:                                    ciliated olfactory receptor neuron
#>  4:                                microvillous olfactory receptor neuron
#>  5:                                       crypt olfactory receptor neuron
#>  6:                                  olfactory epithelial supporting cell
#>  7:                                             olfactory epithelial cell
#>  8:                                    basal cell of olfactory epithelium
#>  9:                                  globose cell of olfactory epithelium
#> 10:                             basal proper cell of olfactory epithelium
#> 11:                                            olfactory ensheathing cell
#> 12:                                            olfactory bulb interneuron
#> 13:                                            olfactory bulb tufted cell
#> 14: nucleus accumbens shell and olfactory tubercle D1 medium spiny neuron
#> 15: nucleus accumbens shell and olfactory tubercle D2 medium spiny neuron
#> 16:                                     tuft cell of olfactory epithelium
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  definition
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      <char>
#>  1:                                                                                                                                                                                                                                                                                                                                                                                                              Any neuron that is capable of some detection of chemical stimulus involved in sensory perception of smell.
#>  2:                                                                                                                                                                   A granule cell that has a soma located in an olfactory bulb granule cell layer. An olfactory granule cell is an interneuron that lacks an axon, makes reciprocal dendro-dendritic synapses with mitral cells and tufted cells and is involved in the fine spatio-temporal tuning of the responses of these principal olfactory bulb neurons to odors.
#>  3:                                                                                                                                                                                                                                                                                                                                                                                                          An olfactory receptor cell in which the apical ending of the dendrite is a pronounced ciliated olfactory knob.
#>  4:                                                                                                                                                                                                                                                                                                                                                                                                         An olfactory receptor cell in which the apical ending of the dendrite is a knob that bears numerous microvilli.
#>  5:                                                                                                                                                                                                                                                                                                                                                                                                                          An olfactory receptor cell with short cilia growing in an invagination bordered by microvilli.
#>  6:               Olfactory epithelial support cell is a columnar cell that extends from the epithelial free margin to the basement membrane of the olfactory epithelium. This cell type has a large, vertically, elongate, euchromatic nucleus, along with other nuclei, forms a layer superficial to the cell body of the receptor cell; sends long somewhat irregular microvilli into the mucus layer; at the base, with expanded end-feet containing numerous lamellated dense bodies resembling lipofuscin of neurons.
#>  7:                                                                                                                                                                                                                                                                                                                                                                                                                                                             A specialized cell involved in sensory perception of smell.
#>  8:                                                                                                                                                                                                                                                                                                                                                                                                                                             An epithelial cell located on the basal lamina of the olfactory epithelium.
#>  9:                                                                                                                                                                                                                                                                                A rounded or elliptical epithelial cell, with pale-staining open face nucleus and pale cytoplasm rich in free ribosomes and clusters of centrioles; form a distinct basal zone spaced slightly from the basal surface of the epithelium.
#> 10:                                                                                                                                                                                                                                                        A flat or angular epithelial cell with condensed nuclei and darkly staining cytoplasm containing numerous intermediate filaments inserted into desmosomes contacting surrounding supporting cells; lie in contact with the basal lamina of olfactory epithelium.
#> 11:                                                                                                                                                                                                                                                                                               A neural-crest derived glial cell that supports the growth and survival of primary olfactory neuroons from the neuroepithelium in the nasal cavity to the brain by encasing large bundles of numerous unmyelinated axons.
#> 12:                                                                                                                                                                                                                                                                                                                                                                                                 A neuron residing in the olfactory bulb that serve to process and refine signals arising from olfactory sensory neurons
#> 13: The principal glutaminergic neuron located in the outer third of the external plexiform layer of the olfactory bulb; a single short primary dendrite traverses the outer external plexiform layer and terminates within an olfactory glomerulus in a tuft of branches, where it receives the input from olfactory receptor neuron axon terminals; axons of the tufted cells transfer information to a number of areas in the brain, including the piriform cortex, entorhinal cortex, olfactory tubercle, and amygdala.
#> 14:                                                                                                                                                                                                                                                                                                                                                                                                                  A DRD1-expressing medium spiny neuron that is part of a nucleus accumbens shell or olfactory tubercle.
#> 15:                                                                                                                                                                                                                                                                                                                                                                                                                  A DRD2-expressing medium spiny neuron that is part of a nucleus accumbens shell or olfactory tubercle.
#> 16:           A tuft cell that is part of the olfactory epithelium, characterized by a globular body and the expression of neurogranin (Nrgn) in mice. This cell plays a crucial role in allergen recognition and regulating olfactory stem cell proliferation via TRPM5-dependent ATP sensing and cysteinyl leukotriene production. Unlike nasal respiratory tuft cells, it has low to absent expression of taste receptors, including the G protein Gα gustducin, and rarely contacts olfactory sensory neurons directly.
```
