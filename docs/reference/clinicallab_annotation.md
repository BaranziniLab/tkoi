# Clinical Laboratory Measurement Annotations (LOINC)

A structured annotation table of clinical laboratory tests and
biomarkers, derived from the LOINC (Logical Observation Identifiers
Names and Codes) system. This dataset maps LOINC identifiers to
human-readable test names and is used within the tKOI pipeline to
annotate graph nodes representing clinical laboratory features.

## Usage

``` r
clinicallab_annotation
```

## Format

A data frame with multiple rows and 2 columns:

- identifier:

  A character string representing the LOINC code (e.g., "2251-7").

- name:

  A human-readable description of the laboratory test or measurement,
  including analyte, method, specimen type, timing, and scale (e.g.,
  "Estriol:MCnc:Pt:Ser/Plas:Qn").

## Source

LOINC — Logical Observation Identifiers Names and Codes
<https://loinc.org/>

## Details

This dataset supports the integration of clinical laboratory data into
biomedical knowledge graphs by providing interpretive labels for
LOINC-coded test results. These annotations are essential for
associating molecular or phenotypic features with clinical diagnostics
in translational research or precision medicine applications.

The dataset is intended to be joined to graph nodes with `identifier`
values corresponding to LOINC codes.

## See also

[`run_tkoi`](run_tkoi.md),
[`disease_annotation`](disease_annotation.md),
[`compound_annotation`](compound_annotation.md)

## Examples

``` r
data(clinicallab_annotation)
subset(clinicallab_annotation, grepl("Estriol", name))
#>       identifier
#> 6         2251-7
#> 1045     27225-2
#> 2253      2250-9
#> 6285      2253-3
#> 6384      2247-5
#> 6888     21265-4
#> 9456     14716-5
#> 14452    34293-1
#> 15956     2249-1
#> 17320    43803-6
#> 19250    47222-5
#> 19441    32353-5
#> 19510    13737-2
#> 21238    60549-3
#> 22046     2252-5
#> 24989    34292-3
#> 27245    44004-0
#> 27275    27259-1
#> 28197    15064-9
#> 32782    20466-9
#> 36461     2248-3
#> 41072    14717-3
#> 44051    43802-8
#> 47276    14718-1
#> 50385    15063-1
#> 53031    16838-5
#> 55152    30510-2
#> 56370    21264-7
#>                                                                                                                                  name
#> 6                                                                                                         Estriol:MCnc:Pt:Ser/Plas:Qn
#> 1045                                                                                                        Estriol:MCnc:Pt:Saliva:Qn
#> 2253                                                                                         Estriol.unconjugated:MCnc:Pt:Ser/Plas:Qn
#> 6285                                                                                                        Estriol:MRat:24H:Urine:Qn
#> 6384                                                                                           Estriol.conjugated:MCnc:Pt:Ser/Plas:Qn
#> 6888                                                                              Estriol.unconjugated^^unadjusted:MoM:Pt:Ser/Plas:Qn
#> 9456                                                                                                      Estriol:SCnc:Pt:Ser/Plas:Qn
#> 14452                                                                                             Estriol/Creatinine:SRto:Pt:Urine:Qn
#> 15956                                                                                       Estriol.unconjugated:MCnc:Pt:Amnio fld:Qn
#> 17320                                                                                            Estriol^^adjusted:MoM:Pt:Ser/Plas:Qn
#> 19250 Trisomy 21 risk based on maternal age + Alpha-1-Fetoprotein + Choriogonadotropin + Estriol.unconjugated:Likelihood:Pt:^Fetus:Qn
#> 19441                                                                                                      Estriol:MoM:Pt:Ser/Plas:Qn
#> 19510                                                                                            Estriol/Creatinine:MRto:24H:Urine:Qn
#> 21238                                                                                                        8-Estriol (U) [Mass/Vol]
#> 22046                                                                                                        Estriol:MCnc:Pt:Urine:Qn
#> 24989                                                                                                       Estriol:SCnc:24H:Urine:Qn
#> 27245                                                                                        Estriol.unconjugated:MoM:Pt:Amnio fld:Qn
#> 27275                                                                                      Estriol.unconjugated:PrThr:Pt:Ser/Plas:Ord
#> 28197                                                                                        Estriol.unconjugated:SCnc:Pt:Ser/Plas:Qn
#> 32782                                                                                         Estriol.unconjugated:MoM:Pt:Ser/Plas:Qn
#> 36461                                                                                                    Estriol:MCnc:Pt:Amnio fld:Qn
#> 41072                                                                                                        Estriol:SCnc:Pt:Urine:Qn
#> 44051                                                                                                       Estriol:MCnc:24H:Urine:Qn
#> 47276                                                                                                       Estriol:SRat:24H:Urine:Qn
#> 50385                                                                                          Estriol.conjugated:SCnc:Pt:Ser/Plas:Qn
#> 53031                                                                                                     Estriol:MCnc:Pt:Body fld:Qn
#> 55152                                                                                             Estriol/Creatinine:MRto:Pt:Urine:Qn
#> 56370                                                                               Estriol.unconjugated^^adjusted:MoM:Pt:Ser/Plas:Qn
```
