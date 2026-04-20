# Protein and Molecular Complex Annotations

A detailed annotation table of known molecular and protein complexes,
including host-pathogen receptor-ligand assemblies and endogenous
regulatory complexes. Each entry provides a unique identifier, a
descriptive name of the complex, and a rich textual definition outlining
biological roles, composition, and functional significance.

## Usage

``` r
complex_annotation
```

## Format

A data frame with 3318 rows and 3 columns:

- identifier:

  A numeric or character string serving as the unique ID for the complex
  (e.g., "2241498559").

- name:

  A descriptive name of the complex, including viral-host receptor
  assemblies or protein transcription factor complexes (e.g.,
  "SARS-CoV-2 Spike - human ACE2 receptor complex").

- description:

  A detailed functional and mechanistic summary of the complex, often
  including molecular interactions, roles in signaling, and implications
  in infection or transcriptional regulation.

## Source

Complex descriptions adapted from curated biological databases such as
Reactome and literature sources.

## Details

This dataset is particularly useful for annotating graph nodes that
represent biological complexes in systems-level analyses. Examples
include virus-receptor binding events relevant to infectious disease
research and SMAD family transcriptional complexes involved in TGF-beta
signaling.

The annotations can be joined to analysis outputs from
[`run_tkoi`](run_tkoi.md) using the `identifier` column.

## See also

[`compound_annotation`](compound_annotation.md),
[`disease_annotation`](disease_annotation.md), [`run_tkoi`](run_tkoi.md)

## Examples

``` r
data(complex_annotation)
subset(complex_annotation, grepl("Spike", name))
#>      identifier                                           name
#> 1    3814064568   SARS-CoV Spike - human ACE2 receptor complex
#> 2    4062616558   SARS-CoV Spike - human CLEC4M lectin complex
#> 3    2241498559 SARS-CoV-2 Spike - human ACE2 receptor complex
#> 4    3898566408  SARS-CoV-2 Spike - human ACE2-SLC6A19 complex
#> 5    2498369513 SARS-CoV-2 Spike - human CLEC4M lectin complex
#> 3318  191534684   MERS-CoV Spike - human DPP4 receptor complex
#>                                                                                                                                                                                                                                                            description
#> 1                                                                                                                                                        Binding of SARS-CoV coronavirus Spike protein to human receptor ACE2 facilitates virus entry into host cell. 
#> 2      Binding of SARS-CoV coronavirus Spike protein to human lectin CLEC4M expressed on innate immune cells and subsequent internalisation may lead to virus clearance or, on the contrary, result in spread of the virus to susceptible cells, or even other hosts. 
#> 3                                                                                                                                                      Binding of SARS-CoV-2 coronavirus Spike protein to human receptor ACE2 facilitates virus entry into host cell. 
#> 4                                                          Internalised SARS-CoV-2 coronavirus Spike - human receptor ACE2 complex binds amino acid transporter SLC6A19 (B0AT1). SLC6A19 may play a regulatory role for the enteric infections of some coronaviruses. 
#> 5    Binding of SARS-CoV-2 coronavirus Spike protein to human lectin CLEC4M expressed on innate immune cells and subsequent internalisation may lead to virus clearance or, on the contrary, result in spread of the virus to susceptible cells, or even other hosts. 
#> 3318                                                                                                                                                         Binding of MERS coronavirus Spike protein to human receptor DPP2 facilitates virus entry into host cell. 
```
