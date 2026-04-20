# Enzyme Commission (EC) Number Annotations

A reference annotation table mapping Enzyme Commission (EC) numbers to
enzyme names. This dataset is used in the tKOI framework to annotate
biological network nodes associated with enzymatic activity.

## Usage

``` r
ec_annotation
```

## Format

A data frame with 8764 rows and 2 columns:

- identifier:

  A character string representing the EC number (e.g., "1.1.1.196"). EC
  numbers classify enzymes by the chemical reactions they catalyze.

- name:

  The full enzyme name (e.g., "15-hydroxyprostaglandin-D dehydrogenase
  (NADP+)"). If missing (`NA`), the enzyme has an identifier but no
  corresponding name in the source.

## Source

IUBMB Enzyme Nomenclature <https://enzyme.expasy.org/>

## Details

EC numbers are a standardized nomenclature for enzymes, maintained by
the IUBMB. This dataset supports functional interpretation of nodes in
the tKOI network with enzymatic roles. The `identifier` column is used
to link to graph nodes of type "EC", and the `name` column provides
descriptive labels for downstream analysis and visualization.

## See also

[`complex_annotation`](complex_annotation.md),
[`run_tkoi`](run_tkoi.md), [`protein_annotation`](protein_annotation.md)

## Examples

``` r
data(ec_annotation)
subset(ec_annotation, grepl("hydroxylase", name, ignore.case = TRUE))
#>       identifier
#> 2     1.14.14.85
#> 8     1.14.15.16
#> 124   1.14.14.71
#> 182    1.17.99.3
#> 189   1.14.13.97
#> 203   1.17.99.10
#> 210   1.14.14.26
#> 313  1.14.14.165
#> 346  1.14.13.172
#> 388   1.14.13.74
#> 482  1.14.14.112
#> 493  1.14.13.159
#> 584   1.14.11.56
#> 733   1.14.14.78
#> 806   1.14.13.73
#> 921   1.14.14.48
#> 933  1.14.14.168
#> 951   1.14.13.99
#> 989  1.14.13.180
#> 1007  1.14.11.39
#> 1146  1.14.11.55
#> 1213  1.14.11.70
#> 1252 1.14.14.146
#> 1253   1.17.99.5
#> 1271  1.14.11.41
#> 1302 1.14.13.191
#> 1312 1.14.14.106
#> 1435  1.14.14.27
#> 1475 1.14.13.240
#> 1502   1.14.99.2
#> 1612 1.14.13.150
#> 1763  1.14.11.40
#> 1783 1.14.14.124
#> 1832  1.14.14.83
#> 1846 1.14.14.145
#> 1854 1.14.14.113
#> 1871  1.14.99.37
#> 1907 1.14.13.152
#> 1933 1.14.13.198
#> 1941  1.14.14.49
#> 1947 1.14.14.178
#> 2109 1.14.14.164
#> 2229 1.14.13.135
#> 2273  1.14.14.29
#> 2360 1.14.13.145
#> 2382  1.14.11.58
#> 2394  1.14.13.94
#> 2404      1.99.1
#> 2415 1.14.13.228
#> 2488  1.14.14.70
#> 2498 1.14.14.158
#> 2512 1.14.14.118
#> 2524 1.14.13.147
#> 2607 1.14.13.220
#> 2634  1.14.11.45
#> 2866  1.14.13.75
#> 2960 1.14.14.121
#> 3000  1.14.11.42
#> 3054  1.14.15.22
#> 3059    1.17.8.1
#> 3063 1.14.13.146
#> 3071 1.14.14.144
#> 3089  1.14.13.63
#> 3176  1.14.13.91
#> 3276  1.14.11.60
#> 3498  1.14.15.27
#> 3514 1.14.14.183
#> 3516    1.17.9.2
#> 3614 1.14.13.219
#> 3630 1.14.14.134
#> 3631 1.14.13.197
#> 3686 1.14.13.181
#> 3729  1.14.14.57
#> 3736 1.14.14.120
#> 3785   1.14.18.8
#> 3845  1.14.11.28
#> 3878  1.14.14.25
#> 3895   1.17.99.2
#> 3931 1.14.13.113
#> 4008  1.14.11.57
#> 4085  1.14.14.54
#> 4090  1.14.14.24
#> 4247 1.14.13.143
#> 4249 1.14.14.182
#> 4359 1.14.13.188
#> 4453  1.14.15.35
#> 4521 1.14.13.119
#> 4588  1.14.99.43
#> 4620  1.14.99.69
#> 4964 1.14.13.123
#> 5020  1.14.11.20
#> 5035  1.14.15.31
#> 5157 1.14.13.100
#> 5176  1.14.11.47
#> 5217  1.14.13.98
#> 5232 1.14.14.136
#> 5262 1.14.13.127
#> 5333  1.14.99.60
#> 5377 1.14.14.139
#> 5395 1.14.14.138
#> 5554 1.14.13.189
#> 5583 1.14.13.109
#> 5752 1.14.13.154
#> 5807 1.14.13.161
#> 5808  1.14.13.77
#> 5895    1.17.3.4
#> 5909  1.14.15.10
#> 6045 1.14.13.209
#> 6054  1.14.15.19
#> 6172    1.17.1.6
#> 6256  1.14.11.61
#> 6278 1.14.14.125
#> 6316  1.14.14.79
#> 6357 1.14.13.199
#> 6409  1.14.99.65
#> 6579   1.17.98.1
#> 6683  1.14.14.95
#> 6745 1.14.13.176
#> 6786  1.14.11.76
#> 6801 1.14.13.153
#> 6809  1.14.15.14
#> 6828  1.14.13.64
#> 6861 1.14.14.103
#> 6961  1.14.11.73
#> 7037 1.14.14.176
#> 7175  1.14.13.95
#> 7239 1.14.13.126
#> 7254 1.14.13.108
#> 7347  1.14.11.35
#> 7353 1.14.14.105
#> 7370  1.14.13.76
#> 7474    1.17.2.2
#> 7511 1.14.14.149
#> 7577     1.7.3.3
#> 7581  1.14.11.71
#> 7677 1.14.14.162
#> 7704 1.14.13.184
#> 7723  1.14.99.49
#> 7736  1.14.14.76
#> 7775 1.14.14.171
#> 7960 1.14.13.194
#> 7970  1.14.11.79
#> 8035  1.14.13.96
#> 8205 1.14.13.110
#> 8233 1.14.14.177
#> 8293 1.14.13.183
#> 8327 1.14.14.104
#> 8364 1.14.14.166
#> 8421 1.14.14.167
#> 8430 1.14.14.163
#> 8449    1.5.99.4
#>                                                                        name
#> 2                                             7-deoxyloganate 7-hydroxylase
#> 8                                                 vitamin D3 24-hydroxylase
#> 124                                          cucurbitadienol 11-hydroxylase
#> 182  3alpha,7alpha,12alpha-trihydroxy-5beta-cholestanoyl-CoA 24-hydroxylase
#> 189                               taurochenodeoxycholate 6alpha-hydroxylase
#> 203                                                steroid C-25 hydroxylase
#> 210                                24-hydroxycholesterol 7alpha-hydroxylase
#> 313                                 indole-3-carbonyl nitrile 4-hydroxylase
#> 346                                                salicylate 5-hydroxylase
#> 388                                            7-deoxyloganin 7-hydroxylase
#> 482                                    ent-cassa-12,15-diene 11-hydroxylase
#> 493                                                vitamin D 25-hydroxylase
#> 584                                             L-proline cis-4-hydroxylase
#> 733                                         phylloquinone omega-hydroxylase
#> 806                                              tabersonine 16-hydroxylase
#> 921                                   jasmonoyl-L-amino acid 12-hydroxylase
#> 933                                     germacrene A acid 8beta-hydroxylase
#> 951                                24-hydroxycholesterol 7alpha-hydroxylase
#> 989                                               aklavinone 12-hydroxylase
#> 1007                                               L-asparagine hydroxylase
#> 1146                                                    ectoine hydroxylase
#> 1213                                  7-deoxycylindrospermopsin hydroxylase
#> 1252                                         geranylgeraniol 18-hydroxylase
#> 1253                                         bile-acid 7alpha-dehydroxylase
#> 1271                                                 L-arginine hydroxylase
#> 1302                                 ent-sandaracopimaradiene 3-hydroxylase
#> 1312                                             taxane 13alpha-hydroxylase
#> 1435                                       resorcinol 4-hydroxylase (FADH2)
#> 1475                                       2-polyprenylphenol 6-hydroxylase
#> 1502                                             kynurenine 7,8-hydroxylase
#> 1612                                          alpha-humulene 10-hydroxylase
#> 1763                                        enduracididine beta-hydroxylase
#> 1783                                         dihydromonacolin L hydroxylase
#> 1832                                                 geraniol 8-hydroxylase
#> 1846                                     abieta-7,13-dien-18-ol hydroxylase
#> 1854                                          alpha-humulene 10-hydroxylase
#> 1871                                           taxadiene 5alpha-hydroxylase
#> 1907                                                 geraniol 8-hydroxylase
#> 1933                                                monacolin L hydroxylase
#> 1941                        12-hydroxyjasmonoyl-L-amino acid 12-hydroxylase
#> 1947                                                steroid 22S-hydroxylase
#> 2109                                                 fraxetin 5-hydroxylase
#> 2229                                     1-hydroxy-2-naphthoate hydroxylase
#> 2273                            25/26-hydroxycholesterol 7alpha-hydroxylase
#> 2360                                   ent-cassa-12,15-diene 11-hydroxylase
#> 2382                        ornithine lipid ester-linked acyl 2-hydroxylase
#> 2394                                         lithocholate 6beta-hydroxylase
#> 2404                                  Hydroxylases (now covered by EC 1.14)
#> 2415                                           jasmonic acid 12-hydroxylase
#> 2488                                 ent-sandaracopimaradiene 3-hydroxylase
#> 2498                                         carotenoid epsilon hydroxylase
#> 2512                                           tryprostatin B 6-hydroxylase
#> 2524                                               taxoid 7beta-hydroxylase
#> 2607                                        resorcinol 4-hydroxylase (NADH)
#> 2634                                             L-isoleucine 4-hydroxylase
#> 2866                                                   vinorine hydroxylase
#> 2960                                          protopanaxadiol 6-hydroxylase
#> 3000         tRNAPhe (7-(3-amino-3-carboxypropyl)wyosine37-C2)-hydroxylase 
#> 3054                                             vitamin D 1,25-hydroxylase
#> 3059                                          hydroxysqualene dehydroxylase
#> 3063                                              taxoid 14beta-hydroxylase
#> 3071                                          abieta-7,13-diene hydroxylase
#> 3089                                   3-hydroxyphenylacetate 6-hydroxylase
#> 3176                                             deoxysarpagine hydroxylase
#> 3276                                               scopoletin 8-hydroxylase
#> 3498                            beta-dihydromenaquinone-9 omega-hydroxylase
#> 3514                                              taxoid 2alpha-hydroxylase
#> 3516                                            (+)-pinoresinol hydroxylase
#> 3614                                       resorcinol 4-hydroxylase (NADPH)
#> 3630                                             beta-amyrin 24-hydroxylase
#> 3631                                         dihydromonacolin L hydroxylase
#> 3686                                       13-deoxydaunorubicin hydroxylase
#> 3729                              taurochenodeoxycholate 6alpha-hydroxylase
#> 3736                                           dammarenediol 12-hydroxylase
#> 3785                   7alpha-hydroxycholest-4-en-3-one 12alpha-hydroxylase
#> 3845                                                  proline 3-hydroxylase
#> 3878                                             cholesterol 24-hydroxylase
#> 3895                                               ethylbenzene hydroxylase
#> 3931                                        FAD-dependent urate hydroxylase
#> 4008                                          L-proline trans-4-hydroxylase
#> 4085                                            phenylacetate 2-hydroxylase
#> 4090                                               vitamin D 25-hydroxylase
#> 4247                                        ent-isokaurene C2/3-hydroxylase
#> 4249                                               taxoid 7beta-hydroxylase
#> 4359                                     6-deoxyerythronolide B hydroxylase
#> 4453                                     6-deoxyerythronolide B hydroxylase
#> 4521                                   5-epiaristolochene 1,3-dihydroxylase
#> 4588                                             beta-amyrin 24-hydroxylase
#> 4620          tRNA 2-(methylsulfanyl)-N6-isopentenyladenosine37 hydroxylase
#> 4964                                               germacrene A hydroxylase
#> 5020                                       deacetoxyvindoline 4-hydroxylase
#> 5035                          2-hydroxy-5-methyl-1-naphthoate 7-hydroxylase
#> 5157                            25/26-hydroxycholesterol 7alpha-hydroxylase
#> 5176                     [50S ribosomal protein L16]-arginine 3-hydroxylase
#> 5217                                             cholesterol 24-hydroxylase
#> 5232                                             deoxysarpagine hydroxylase
#> 5262                              3-(3-hydroxyphenyl)propanoate hydroxylase
#> 5333                                     3-demethoxyubiquinol 3-hydroxylase
#> 5377                5beta-cholestane-3alpha,7alpha-diol 12alpha-hydroxylase
#> 5395                                         lithocholate 6beta-hydroxylase
#> 5554                                    5-methyl-1-naphthoate 3-hydroxylase
#> 5583                                     abieta-7,13-dien-18-ol hydroxylase
#> 5752                                            erythromycin 12-hydroxylase
#> 5807                                          (+)-camphor 6-exo-hydroxylase
#> 5808                                             taxane 13alpha-hydroxylase
#> 5895                                                  juglone 3-hydroxylase
#> 5909                                         (+)-camphor 6-endo-hydroxylase
#> 6045                                           salicyloyl-CoA 5-hydroxylase
#> 6054                                        C-19 steroid 1alpha-hydroxylase
#> 6172                                         bile-acid 7alpha-dehydroxylase
#> 6256                                             feruloyl-CoA 6-hydroxylase
#> 6278                                                monacolin L hydroxylase
#> 6316                                 docosahexaenoic acid omega-hydroxylase
#> 6357                                 docosahexaenoic acid omega-hydroxylase
#> 6409   4-amino-L-phenylalanyl-[CmlP-peptidyl-carrier-protein] 3-hydroxylase
#> 6579                                         bile-acid 7alpha-dehydroxylase
#> 6683                                               germacrene A hydroxylase
#> 6745                                           tryprostatin B 6-hydroxylase
#> 6786                                           L-glutamate 3(R)-hydroxylase
#> 6801                                             (+)-sabinene 3-hydroxylase
#> 6809                                methyl-branched lipid omega-hydroxylase
#> 6828                                        4-hydroxybenzoate 1-hydroxylase
#> 6861                                             tabersonine 16-hydroxylase
#> 6961                                       [protein]-arginine 3-hydroxylase
#> 7037                                           taxadiene 5alpha-hydroxylase
#> 7175                   7alpha-hydroxycholest-4-en-3-one 12alpha-hydroxylase
#> 7239                                              vitamin D3 24-hydroxylase
#> 7254                                          abieta-7,13-diene hydroxylase
#> 7347                              1-deoxypentalenic acid 11beta-hydroxylase
#> 7353                                              taxane 10beta-hydroxylase
#> 7370                                              taxane 10beta-hydroxylase
#> 7474                                 lupanine 17-hydroxylase (cytochrome c)
#> 7511                                   5-epiaristolochene 1,3-dihydroxylase
#> 7577                                   factor-independent urate hydroxylase
#> 7581                                          methylphosphonate hydroxylase
#> 7677                                                flavanone 2-hydroxylase
#> 7704                                          protopanaxadiol 6-hydroxylase
#> 7723                          2-hydroxy-5-methyl-1-naphthoate 7-hydroxylase
#> 7736                                       ent-isokaurene C2/C3-hydroxylase
#> 7775                                        beta-amyrin 16alpha-hydroxylase
#> 7960                                        phylloquinone omega-hydroxylase
#> 7970                                 protein-L-histidine (3S)-3-hydroxylase
#> 8035                5beta-cholestane-3alpha,7alpha-diol 12alpha-hydroxylase
#> 8205                                         geranylgeraniol 18-hydroxylase
#> 8233                          ultra-long-chain fatty acid omega-hydroxylase
#> 8293                                           dammarenediol 12-hydroxylase
#> 8327                                                   vinorine hydroxylase
#> 8364                                     (S)-N-methylcanadine 1-hydroxylase
#> 8421         (13S,14R)-13-O-acetyl-1-hydroxy-N-methylcanadine 8-hydroxylase
#> 8430                          (S)-1-hydroxy-N-methylcanadine 13-hydroxylase
#> 8449                                                 nicotine 6-hydroxylase
```
