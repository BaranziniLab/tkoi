# Protein Domain Annotations (Pfam)

A curated annotation table of protein domains based on Pfam identifiers.
This dataset provides structured mappings from Pfam domain IDs to domain
names, supporting the annotation and interpretation of domain-level
biological features in knowledge graphs used by the tKOI framework.

## Usage

``` r
proteindomain_annotation
```

## Format

A data frame with 14193 rows and 2 columns:

- identifier:

  A character string representing the Pfam domain identifier (e.g.,
  "PF19543").

- name:

  The name or description of the domain, often including structural or
  functional attributes (e.g., "Glycoside hydrolase 123, N-terminal
  domain").

## Source

Pfam Database <https://pfam.xfam.org/>

## Details

Protein domains represent conserved functional or structural units
within proteins. This dataset enables enrichment and network-based
inference at the domain level. The `identifier` column can be joined to
graph nodes representing domains in tKOI output, while the `name` column
adds interpretability for visualization and reporting.

Domains include a variety of structural motifs, enzyme catalytic
regions, and uncharacterized conserved segments (e.g., DUFs).

## See also

[`protein_annotation`](protein_annotation.md),
[`proteinfamily_annotation`](proteinfamily_annotation.md),
[`run_tkoi`](run_tkoi.md)

## Examples

``` r
data(proteindomain_annotation)
subset(proteindomain_annotation, grepl("hydrolase", name, ignore.case = TRUE))
#>       identifier
#> 1        PF19543
#> 46       PF00657
#> 144      PF00933
#> 170      PF08307
#> 224      PF02964
#> 231      PF02275
#> 356      PF13318
#> 362      PF01670
#> 411      PF03512
#> 470      PF03636
#> 610      PF14508
#> 648      PF01915
#> 667      PF14323
#> 788      PF02638
#> 793      PF02550
#> 827      PF09370
#> 897      PF19420
#> 930      PF16923
#> 978      PF06342
#> 980      PF00443
#> 1101     PF16317
#> 1226     PF18271
#> 1278     PF00670
#> 1303     PF01227
#> 1311     PF13647
#> 1322     PF12471
#> 1343     PF02156
#> 1352     PF10605
#> 1417     PF16255
#> 1519     PF14196
#> 1583     PF13869
#> 1584     PF03747
#> 1667     PF07745
#> 1689     PF04273
#> 1700     PF18089
#> 1928     PF00723
#> 1931     PF05908
#> 1993     PF00722
#> 2006     PF17433
#> 2047     PF02011
#> 2114     PF16875
#> 2131     PF12215
#> 2190     PF02057
#> 2247     PF14587
#> 2340     PF03662
#> 2378     PF01183
#> 2491     PF00332
#> 2532     PF04996
#> 2552     PF06439
#> 2582     PF09752
#> 2618     PF20091
#> 2619     PF14509
#> 2647     PF00704
#> 2673     PF13423
#> 2868     PF03663
#> 3098     PF17167
#> 3153     PF02633
#> 3171     PF05838
#> 3185     PF01074
#> 3480     PF12697
#> 3509     PF12535
#> 3515     PF11308
#> 3652     PF08282
#> 3658     PF18565
#> 3832     PF18034
#> 3922     PF02015
#> 4000     PF08472
#> 4046     PF10566
#> 4078     PF10340
#> 4143     PF15979
#> 4146     PF07858
#> 4157     PF16862
#> 4172     PF18564
#> 4178     PF03537
#> 4193     PF02882
#> 4302     PF11790
#> 4314     PF01156
#> 4320     PF07461
#> 4326     PF00925
#> 4327     PF05165
#> 4429     PF07488
#> 4459     PF03403
#> 4524     PF02148
#> 4547     PF01503
#> 4610     PF01981
#> 4658     PF12891
#> 4660     PF07826
#> 4781     PF07859
#> 4865     PF14885
#> 4917     PF02837
#> 4963     PF17189
#> 4994     PF08306
#> 5028     PF14871
#> 5126     PF19718
#> 5328     PF02838
#> 5366     PF00150
#> 5390     PF01738
#> 5490     PF08244
#> 5628     PF17451
#> 5705     PF14883
#> 5801     PF03648
#> 6030     PF13419
#> 6102     PF13199
#> 6130     PF06259
#> 6132     PF09081
#> 6141     PF00728
#> 6154     PF03959
#> 6164     PF02836
#> 6204     PF06399
#> 6233     PF12695
#> 6237     PF00702
#> 6754     PF16123
#> 6755     PF17387
#> 6763     PF03664
#> 6802     PF04685
#> 6880     PF01301
#> 6967     PF01270
#> 6981     PF12715
#> 7028     PF00232
#> 7165     PF01374
#> 7255     PF03659
#> 7278     PF18031
#> 7435     PF12888
#> 7500     PF14736
#> 7576     PF04909
#> 7620     PF00457
#> 7686     PF17678
#> 7694     PF03632
#> 7697     PF13336
#> 7756     PF01088
#> 7948     PF02649
#> 8047     PF07521
#> 8106     PF02055
#> 8178     PF04616
#> 8189     PF17890
#> 8397     PF07477
#> 8484     PF01195
#> 8561     PF02324
#> 8636     PF08840
#> 8652     PF17652
#> 8694     PF12710
#> 8978     PF01979
#> 9013     PF09127
#> 9050     PF17677
#> 9124     PF03065
#> 9176     PF07486
#> 9192     PF05116
#> 9217     PF10230
#> 9506     PF03633
#> 9716     PF14872
#> 9724     PF03200
#> 9729     PF18230
#> 9732     PF03644
#> 9766     PF00295
#> 9805     PF18088
#> 9980     PF11975
#> 9991     PF03639
#> 10035    PF16874
#> 10057    PF16674
#> 10089    PF02289
#> 10090    PF03718
#> 10203    PF13472
#> 10273    PF07969
#> 10336    PF01557
#> 10408    PF01229
#> 10445    PF00703
#> 10481    PF05221
#> 10559    PF05013
#> 10619    PF00251
#> 10626    PF07176
#> 10659    PF20408
#> 10842    PF07335
#> 10878    PF10081
#> 11021    PF16822
#> 11038    PF09663
#> 11073    PF00795
#> 11088    PF00331
#> 11215    PF00561
#> 11263    PF09994
#> 11385    PF13872
#> 11420    PF03819
#> 11586    PF01341
#> 11683    PF00840
#> 11693    PF01532
#> 11700    PF01373
#> 11838    PF18438
#> 11913    PF00759
#> 11916    PF02056
#> 11935    PF13286
#> 12016    PF18290
#> 12044    PF14498
#> 12103    PF04447
#> 12158    PF13344
#> 12201    PF07748
#> 12302    PF04083
#> 12343    PF20309
#> 12352    PF04307
#> 12441    PF06441
#> 12459    PF17829
#> 12516    PF06821
#> 12597    PF07698
#> 12653    PF00763
#> 12724    PF01502
#> 12727    PF13200
#> 12854    PF04775
#> 13003    PF05990
#> 13019    PF07470
#> 13065    PF06028
#> 13075    PF14606
#> 13085    PF05382
#> 13161    PF07971
#> 13247    PF21365
#> 13258    PF07632
#> 13270    PF21246
#> 13290    PF08924
#> 13388    PF13320
#> 13426    PF16141
#> 13471    PF13802
#> 13519    PF21252
#> 13526    PF20811
#> 13529    PF21570
#> 13651    PF21087
#> 13658    PF05028
#> 13875    PF01055
#> 13941    PF21104
#> 14018    PF16347
#>                                                                       name
#> 1                               Glycoside hydrolase 123, N-terminal domain
#> 46                                          GDSL-like Lipase/Acylhydrolase
#> 144                          Glycosyl hydrolase family 3 N terminal domain
#> 170                         Glycosyl hydrolase family 98 C-terminal domain
#> 224                           Methane monooxygenase, hydrolase gamma chain
#> 231           Linear amide C-N hydrolases, choloylglycine hydrolase family
#> 356                            1-carboxybiuret hydrolase subunit AtzG-like
#> 362                                           Glycosyl hydrolase family 12
#> 411                                           Glycosyl hydrolase family 52
#> 470                       Glycosyl hydrolase family 65, N-terminal domain 
#> 610                                       Glycosyl-hydrolase 97 N-terminal
#> 648                          Glycosyl hydrolase family 3 C-terminal domain
#> 667                 GxGYxYP putative glycoside hydrolase C-terminal domain
#> 788                                             Glycosyl hydrolase-like 10
#> 793                     Acetyl-CoA hydrolase/transferase N-terminal domain
#> 827                                     Phosphoenolpyruvate hydrolase-like
#> 897                     N,N dimethylarginine dimethylhydrolase, eukaryotic
#> 930                         Glycosyl hydrolase family 63 N-terminal domain
#> 978                     Alpha/beta hydrolase of unknown function (DUF1057)
#> 980                                  Ubiquitin carboxyl-terminal hydrolase
#> 1101                                          Glycosyl hydrolase family 99
#> 1226                   Glycoside hydrolase 131 catalytic N-terminal domain
#> 1278               S-adenosyl-L-homocysteine hydrolase, NAD binding domain
#> 1303                                                  GTP cyclohydrolase I
#> 1311                         Glycosyl hydrolase family 80 of chitosanase A
#> 1322                                        GTP cyclohydrolase N terminal 
#> 1343                                          Glycosyl hydrolase family 26
#> 1352                                       3HB-oligomer hydrolase (3HBOH) 
#> 1417                                        GDSL-like Lipase/Acylhydrolase
#> 1519                      L-2-amino-thiazoline-4-carboxylic acid hydrolase
#> 1583                                                  Nucleotide hydrolase
#> 1584                                             ADP-ribosylglycohydrolase
#> 1667                                          Glycosyl hydrolase family 53
#> 1689        Beta-lactamase hydrolase-like protein, phosphatase-like domain
#> 1700                                            DAPG hydrolase PhiG domain
#> 1928                                         Glycosyl hydrolases family 15
#> 1931                                        Poly-gamma-glutamate hydrolase
#> 1993                                         Glycosyl hydrolases family 16
#> 2006                Glycosyl hydrolase family 49 N-terminal Ig-like domain
#> 2047                                          Glycosyl hydrolase family 48
#> 2114                        Glycosyl hydrolase family 36 N-terminal domain
#> 2131              beta-glucosidase 2, glycosyl-hydrolase family 116 N-term
#> 2190                                          Glycosyl hydrolase family 59
#> 2247                                        O-Glycosyl hydrolase family 30
#> 2340                      Glycosyl hydrolase family 79, N-terminal domain 
#> 2378                                         Glycosyl hydrolases family 25
#> 2491                                         Glycosyl hydrolases family 17
#> 2532                                          Succinylarginine dihydrolase
#> 2552                                         3-keto-disaccharide hydrolase
#> 2582                            Alpha/beta hydrolase domain containing 18 
#> 2618                                           Alpha/beta hydrolase domain
#> 2619                     Glycosyl-hydrolase 97 C-terminal, oligomerisation
#> 2647                                         Glycosyl hydrolases family 18
#> 2673                                 Ubiquitin carboxyl-terminal hydrolase
#> 2868                                         Glycosyl hydrolase family 76 
#> 3098                   Glycosyl hydrolase 36 superfamily, catalytic domain
#> 3153                                             Creatinine amidohydrolase
#> 3171                                                Glycosyl hydrolase 108
#> 3185                       Glycosyl hydrolases family 38 N-terminal domain
#> 3480                                           Alpha/beta hydrolase family
#> 3509              Hydrolase of X-linked nucleoside diphosphate N terminal 
#> 3515                    Glycosyl hydrolases related to GH101 family, GH129
#> 3652                                  haloacid dehalogenase-like hydrolase
#> 3658                      Glycoside hydrolase family 2 C-terminal domain 5
#> 3832               Bacterial Glycosyl hydrolase family 3 C-terminal domain
#> 3922                                          Glycosyl hydrolase family 45
#> 4000                       Sucrose-6-phosphate phosphohydrolase C-terminal
#> 4046                                              Glycoside hydrolase 97  
#> 4078                                               Steryl acetyl hydrolase
#> 4143                                         Glycosyl hydrolase family 115
#> 4146                       Limonene-1,2-epoxide hydrolase catalytic domain
#> 4157                   Glycosyl hydrolase family 79 C-terminal beta domain
#> 4172                        Glycoside hydrolase family 5 C-terminal domain
#> 4178                                      Glycoside-hydrolase family GH114
#> 4193  Tetrahydrofolate dehydrogenase/cyclohydrolase, NAD(P)-binding domain
#> 4302                                     Glycosyl hydrolase catalytic core
#> 4314                       Inosine-uridine preferring nucleoside hydrolase
#> 4320                 Nicotine adenine dinucleotide glycohydrolase (NADase)
#> 4326                                                 GTP cyclohydrolase II
#> 4327                                                GTP cyclohydrolase III
#> 4429                            Glycosyl hydrolase family 67 middle domain
#> 4459                Platelet-activating factor acetylhydrolase, isoform II
#> 4524                   Zn-finger in ubiquitin-hydrolases and other protein
#> 4547                               Phosphoribosyl-ATP pyrophosphohydrolase
#> 4610                                          Peptidyl-tRNA hydrolase PTH2
#> 4658                                         Glycoside hydrolase family 44
#> 4660                                       IMP cyclohydrolase-like protein
#> 4781                                             alpha/beta hydrolase fold
#> 4865                             Hypothetical glycosyl hydrolase family 15
#> 4917                    Glycosyl hydrolases family 2, sugar binding domain
#> 4963                     Glycosyl hydrolase family 30 beta sandwich domain
#> 4994                                          Glycosyl hydrolase family 98
#> 5028                                     Hypothetical glycosyl hydrolase 6
#> 5126                   Ubiquitin carboxyl-terminal hydrolase 47 C-terminal
#> 5328                                Glycosyl hydrolase family 20, domain 2
#> 5366                               Cellulase (glycosyl hydrolase family 5)
#> 5390                                         Dienelactone hydrolase family
#> 5490                              Glycosyl hydrolases family 32 C terminal
#> 5628                           Glycosyl hydrolase 101 beta sandwich domain
#> 5705                             Hypothetical glycosyl hydrolase family 13
#> 5801                               Glycosyl hydrolase family 67 N-terminus
#> 6030                                  Haloacid dehalogenase-like hydrolase
#> 6102                                          Glycosyl hydrolase family 66
#> 6130                                                  Alpha/beta hydrolase
#> 6132                       Glucan 1,4-alpha-maltotetraohydrolase, domain C
#> 6141                        Glycosyl hydrolase family 20, catalytic domain
#> 6154                                               Serine hydrolase (FSH1)
#> 6164                       Glycosyl hydrolases family 2, TIM barrel domain
#> 6204               GTP cyclohydrolase I feedback regulatory protein (GFRP)
#> 6233                                           Alpha/beta hydrolase family
#> 6237                                  haloacid dehalogenase-like hydrolase
#> 6754                           Hydroxyacylglutathione hydrolase C-terminus
#> 6755                           Glycosyl hydrolase family 59 central domain
#> 6763                                         Glycosyl hydrolase family 62 
#> 6802                       Glycosyl-hydrolase family 116, catalytic region
#> 6880                                         Glycosyl hydrolases family 35
#> 6967                                          Glycosyl hydrolases family 8
#> 6981                                                    Abhydrolase family
#> 7028                                           Glycosyl hydrolase family 1
#> 7165                                          Glycosyl hydrolase family 46
#> 7255                                          Glycosyl hydrolase family 71
#> 7278                                Ubiquitin carboxyl-terminal hydrolases
#> 7435                                      Lipid-binding putative hydrolase
#> 7500                          Protein N-terminal asparagine amidohydrolase
#> 7576                                                        Amidohydrolase
#> 7620                                         Glycosyl hydrolases family 11
#> 7686                        Glycosyl hydrolase family 92 N-terminal domain
#> 7694                 Glycosyl hydrolase family 65 central catalytic domain
#> 7697                    Acetyl-CoA hydrolase/transferase C-terminal domain
#> 7756                       Ubiquitin carboxyl-terminal hydrolase, family 1
#> 7948                                       Type I GTP cyclohydrolase folE2
#> 8047                 Zn-dependent metallo-hydrolase RNA specificity domain
#> 8106                        Glycosyl hydrolase family 30 TIM-barrel domain
#> 8178                                         Glycosyl hydrolases family 43
#> 8189                          Peptidoglycan hydrolase LytB WW-like domain 
#> 8397                               Glycosyl hydrolase family 67 C-terminus
#> 8484                                               Peptidyl-tRNA hydrolase
#> 8561                                          Glycosyl hydrolase family 70
#> 8636                        BAAT / Acyl-CoA thioester hydrolase C terminal
#> 8652                        Glycosyl hydrolase family 81 C-terminal domain
#> 8694                                  haloacid dehalogenase-like hydrolase
#> 8978                                                 Amidohydrolase family
#> 9013                                  Leukotriene A4 hydrolase, C-terminal
#> 9050         Glycosyl hydrolases family 38 C-terminal beta sandwich domain
#> 9124                                          Glycosyl hydrolase family 57
#> 9176                                                   Cell Wall Hydrolase
#> 9192                                 Sucrose-6F-phosphate phosphohydrolase
#> 9217                                    Lipid-droplet associated hydrolase
#> 9506                      Glycosyl hydrolase family 65, C-terminal domain 
#> 9716                                    Hypothetical glycoside hydrolase 5
#> 9724                        Glycosyl hydrolase family 63 C-terminal domain
#> 9729                   Glycosyl hydrolases family 38 C-terminal sub-domain
#> 9732                                         Glycosyl hydrolase family 85 
#> 9766                                         Glycosyl hydrolases family 28
#> 9805                             Glycoside Hydrolase 20C C-terminal domain
#> 9980                         Family 4 glycosyl hydrolase C-terminal domain
#> 9991                        Glycosyl hydrolase family 81 N-terminal domain
#> 10035                       Glycosyl hydrolase family 36 C-terminal domain
#> 10057               N-terminal of ubiquitin carboxyl-terminal hydrolase 37
#> 10089                                                 Cyclohydrolase (MCH)
#> 10090                                         Glycosyl hydrolase family 49
#> 10203                                GDSL-like Lipase/Acylhydrolase family
#> 10273                                                Amidohydrolase family
#> 10336                           Fumarylacetoacetate (FAA) hydrolase family
#> 10408                                        Glycosyl hydrolases family 39
#> 10445                                         Glycosyl hydrolases family 2
#> 10481                                  S-adenosyl-L-homocysteine hydrolase
#> 10559                                     N-formylglutamate amidohydrolase
#> 10619                      Glycosyl hydrolases family 32 N-terminal domain
#> 10626                   Alpha/beta hydrolase of unknown function (DUF1400)
#> 10659                                          Alpha/beta hydrolase domain
#> 10842                    Fungal chitosanase of glycosyl hydrolase group 75
#> 10878                                          Alpha/beta-hydrolase family
#> 11021                   SGNH hydrolase-like domain, acetyltransferase AlgX
#> 11038                Amidohydrolase ring-opening protein (Amido_AtzD_TrzD)
#> 11073                                            Carbon-nitrogen hydrolase
#> 11088                                         Glycosyl hydrolase family 10
#> 11215                                            alpha/beta hydrolase fold
#> 11263                Uncharacterized alpha/beta hydrolase domain (DUF2235)
#> 11385                               P-loop containing NTP hydrolase pore-1
#> 11420                          MazG nucleotide pyrophosphohydrolase domain
#> 11586                                         Glycosyl hydrolases family 6
#> 11683                                          Glycosyl hydrolase family 7
#> 11693                                         Glycosyl hydrolase family 47
#> 11700                                         Glycosyl hydrolase family 14
#> 11838                    Glycosyl hydrolases family 38 C-terminal domain 1
#> 11913                                          Glycosyl hydrolase family 9
#> 11916                                          Family 4 glycosyl hydrolase
#> 11935                                   Phosphohydrolase-associated domain
#> 12016                                               Nudix hydrolase domain
#> 12044                      Glycosyl hydrolase family 65, N-terminal domain
#> 12103                                       dATP/dGTP pyrophosphohydrolase
#> 12158                                 Haloacid dehalogenase-like hydrolase
#> 12201                      Glycosyl hydrolases family 38 C-terminal domain
#> 12302                           Partial alpha/beta-hydrolase lipase region
#> 12343             Deoxyribohydrolase (DRHyd) domain of the ASK signalosome
#> 12352           LexA-binding, inner membrane-associated putative hydrolase
#> 12441                                         Epoxide hydrolase N terminus
#> 12459                      Gylcosyl hydrolase family 115 C-terminal domain
#> 12516                                                     Serine hydrolase
#> 12597                         7TM receptor with intracellular HD hydrolase
#> 12653      Tetrahydrofolate dehydrogenase/cyclohydrolase, catalytic domain
#> 12724                                    Phosphoribosyl-AMP cyclohydrolase
#> 12727                                   Putative glycosyl hydrolase domain
#> 12854                  Acyl-CoA thioester hydrolase/BAAT N-terminal region
#> 13003                    Alpha/beta hydrolase of unknown function (DUF900)
#> 13019                                         Glycosyl Hydrolase Family 88
#> 13065                    Alpha/beta hydrolase of unknown function (DUF915)
#> 13075                                GDSL-like Lipase/Acylhydrolase family
#> 13085                               Bacteriophage peptidoglycan hydrolase 
#> 13161                        Glycosyl hydrolase family 92 catalytic domain
#> 13247                       Glycosyl hydrolase family 31 C-terminal domain
#> 13258           Cellulose-binding Sde182, nucleoside hydrolase-like domain
#> 13270     Ubiquitin carboxyl-terminal hydrolase 38-like, N-terminal domain
#> 13290                        Rv2525c-like, glycoside hydrolase-like domain
#> 13388                            Glycoside hydrolase 123, catalytic domain
#> 13426                           Glycoside hydrolase family 18, BT1044-like
#> 13471    Glycosyl hydrolase 31 N-terminal galactose mutarotase-like domain
#> 13519                            Glycosyl hydrolase 109, C-terminal domain
#> 13526              Poly (ADP-ribose) glycohydrolase (PARG), helical domain
#> 13529            Arginine dihydrolase ArgZ-like, C-terminal, Rossmann fold
#> 13651                                        Glycosyl hydrolase family 134
#> 13658           Poly (ADP-ribose) glycohydrolase (PARG), Macro domain fold
#> 13875                      Glycosyl hydrolases family 31 TIM-barrel domain
#> 13941    Glycosyl hydrolase family 78 alpha-rhamnosidase N-terminal domain
#> 14018                      N-sulphoglucosamine sulphohydrolase, C-terminal
```
