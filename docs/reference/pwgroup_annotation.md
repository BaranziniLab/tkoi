# Pathway Group Annotations (Reactome Interaction Sets)

A curated annotation table for grouped pathway interaction events from
Reactome. Each entry represents a molecular interaction or signaling
complex involved in a pathway context, typically grouped by protein
complexes or ligand-receptor interactions. These pathway groupings
support detailed annotation of graph nodes in the tKOI network that
represent composite signaling units.

## Usage

``` r
pwgroup_annotation
```

## Format

A data frame with 6343 rows and 2 columns:

- identifier:

  A character string representing the Reactome stable identifier (e.g.,
  "reactome:R-HSA-9624419").

- name:

  A human-readable description of the interaction group or signaling
  unit (e.g., "EGFR:EGF-like ligands").

## Source

Reactome Pathway Database <https://reactome.org>

## Details

This dataset provides molecular resolution for pathway-associated
interactions and groupings (e.g., dimerized receptors, ligand-bound
complexes), enhancing the interpretability of nodes labeled as "PwGroup"
in tKOI analyses. These groups often represent intermediate or
context-specific biochemical states within signaling cascades.

The `identifier` field can be used to join with graph nodes, while the
`name` field is useful for labeling and summarizing pathway-level
results.

## See also

[`pathway_annotation`](pathway_annotation.md),
[`reaction_annotation`](reaction_annotation.md),
[`run_tkoi`](run_tkoi.md)

## Examples

``` r
data(pwgroup_annotation)
subset(pwgroup_annotation, grepl("EGFR", name))
#>                  identifier
#> 2    reactome:R-HSA-9674862
#> 3    reactome:R-HSA-9624419
#> 4    reactome:R-HSA-9624420
#> 5    reactome:R-HSA-9624425
#> 6     reactome:R-HSA-212717
#> 7     reactome:R-HSA-212703
#> 9     reactome:R-HSA-179820
#> 12   reactome:R-HSA-9624427
#> 13    reactome:R-HSA-180337
#> 14    reactome:R-HSA-180331
#> 18    reactome:R-HSA-180348
#> 19    reactome:R-HSA-180286
#> 24    reactome:R-HSA-179791
#> 25    reactome:R-HSA-180269
#> 26    reactome:R-HSA-180503
#> 30    reactome:R-HSA-180326
#> 31    reactome:R-HSA-180343
#> 32    reactome:R-HSA-180277
#> 35    reactome:R-HSA-182960
#> 36    reactome:R-HSA-182953
#> 37    reactome:R-HSA-182939
#> 39    reactome:R-HSA-182930
#> 41    reactome:R-HSA-182928
#> 42    reactome:R-HSA-182948
#> 43    reactome:R-HSA-182945
#> 49    reactome:R-HSA-182946
#> 56   reactome:R-HSA-8951489
#> 57    reactome:R-HSA-182932
#> 60    reactome:R-HSA-182961
#> 61    reactome:R-HSA-182935
#> 62    reactome:R-HSA-182936
#> 65    reactome:R-HSA-182931
#> 66    reactome:R-HSA-179882
#> 67    reactome:R-HSA-179860
#> 68   reactome:R-HSA-8864028
#> 69   reactome:R-HSA-8864018
#> 70   reactome:R-HSA-8867034
#> 71   reactome:R-HSA-8867035
#> 531   reactome:R-HSA-195398
#> 534   reactome:R-HSA-195403
#> 537   reactome:R-HSA-195405
#> 539   reactome:R-HSA-195402
#> 541  reactome:R-HSA-4420101
#> 542  reactome:R-HSA-4420096
#> 543  reactome:R-HSA-4420196
#> 544  reactome:R-HSA-4420164
#> 545  reactome:R-HSA-4420141
#> 546  reactome:R-HSA-5218639
#> 547  reactome:R-HSA-5218641
#> 548  reactome:R-HSA-5218636
#> 549  reactome:R-HSA-5218635
#> 551  reactome:R-HSA-5218779
#> 552  reactome:R-HSA-5218637
#> 557  reactome:R-HSA-5218775
#> 558  reactome:R-HSA-5218771
#> 559  reactome:R-HSA-5218787
#> 560  reactome:R-HSA-5218778
#> 561  reactome:R-HSA-5218757
#> 562  reactome:R-HSA-5218760
#> 563  reactome:R-HSA-5218781
#> 564  reactome:R-HSA-5218756
#> 565  reactome:R-HSA-5218790
#> 566  reactome:R-HSA-5218762
#> 567  reactome:R-HSA-5218759
#> 568  reactome:R-HSA-5218765
#> 570  reactome:R-HSA-5218776
#> 574  reactome:R-HSA-4420171
#> 575  reactome:R-HSA-4420137
#> 576  reactome:R-HSA-4420169
#> 577  reactome:R-HSA-5218788
#> 578  reactome:R-HSA-5357442
#> 579  reactome:R-HSA-5357462
#> 580  reactome:R-HSA-5357453
#> 581  reactome:R-HSA-5218777
#> 615  reactome:R-HSA-5218798
#> 616  reactome:R-HSA-5218786
#> 617  reactome:R-HSA-5218795
#> 618  reactome:R-HSA-5218780
#> 619  reactome:R-HSA-5218766
#> 622  reactome:R-HSA-5218782
#> 624  reactome:R-HSA-5218792
#> 627  reactome:R-HSA-5218793
#> 631  reactome:R-HSA-4420203
#> 632  reactome:R-HSA-4420110
#> 643   reactome:R-HSA-195419
#> 644   reactome:R-HSA-195427
#> 688  reactome:R-HSA-1963571
#> 689   reactome:R-HSA-179847
#> 695  reactome:R-HSA-1227939
#> 701  reactome:R-HSA-1810434
#> 705  reactome:R-HSA-1248703
#> 706  reactome:R-HSA-1248702
#> 713  reactome:R-HSA-1963593
#> 728  reactome:R-HSA-1963587
#> 729  reactome:R-HSA-1250505
#> 736  reactome:R-HSA-1306958
#> 737  reactome:R-HSA-1306961
#> 738  reactome:R-HSA-1306962
#> 739  reactome:R-HSA-1251942
#> 749  reactome:R-HSA-8863840
#> 752  reactome:R-HSA-8863834
#> 755  reactome:R-HSA-8864106
#> 756  reactome:R-HSA-8864118
#> 757  reactome:R-HSA-8864120
#> 758  reactome:R-HSA-8864124
#> 759  reactome:R-HSA-8864121
#> 778  reactome:R-HSA-1977956
#> 1532 reactome:R-HSA-2179410
#> 1533 reactome:R-HSA-2179388
#> 1534 reactome:R-HSA-2179409
#> 1729 reactome:R-HSA-9018569
#> 2256 reactome:R-HSA-9625453
#> 2257 reactome:R-HSA-9625464
#> 2941 reactome:R-HSA-8857547
#> 2942 reactome:R-HSA-8857546
#> 2943 reactome:R-HSA-8857556
#> 2944 reactome:R-HSA-8857567
#> 2945 reactome:R-HSA-8857574
#> 3135 reactome:R-HSA-1218825
#> 3136 reactome:R-HSA-1220581
#> 3137 reactome:R-HSA-1220583
#> 3138 reactome:R-HSA-1220578
#> 3139 reactome:R-HSA-1220573
#> 3140 reactome:R-HSA-1220580
#> 3141 reactome:R-HSA-1500847
#> 3142 reactome:R-HSA-1500849
#> 3143 reactome:R-HSA-1220582
#> 3144 reactome:R-HSA-1181062
#> 3145 reactome:R-HSA-1220652
#> 3146 reactome:R-HSA-1220648
#> 3147 reactome:R-HSA-1809164
#> 3148 reactome:R-HSA-1809167
#> 3149 reactome:R-HSA-1220654
#> 3150 reactome:R-HSA-1220650
#> 3151 reactome:R-HSA-1809168
#> 3152 reactome:R-HSA-1809170
#> 3153 reactome:R-HSA-1225859
#> 3154 reactome:R-HSA-1225858
#> 3155 reactome:R-HSA-1225854
#> 3156 reactome:R-HSA-1225856
#> 3157 reactome:R-HSA-1226015
#> 3158 reactome:R-HSA-1226013
#> 3159 reactome:R-HSA-1247845
#> 3160 reactome:R-HSA-1247846
#> 3161 reactome:R-HSA-1225855
#> 3162 reactome:R-HSA-1225857
#> 3164 reactome:R-HSA-1220585
#> 3165 reactome:R-HSA-1220574
#> 3166 reactome:R-HSA-1220575
#> 3167 reactome:R-HSA-1220587
#> 3168 reactome:R-HSA-1220572
#> 3169 reactome:R-HSA-1220571
#> 3170 reactome:R-HSA-1181059
#> 3171 reactome:R-HSA-1220579
#> 3172 reactome:R-HSA-1248004
#> 3173 reactome:R-HSA-1248010
#> 3174 reactome:R-HSA-1248654
#> 3175 reactome:R-HSA-1248651
#> 3176 reactome:R-HSA-5637697
#> 3177 reactome:R-HSA-5637701
#> 3178 reactome:R-HSA-5637703
#> 3179 reactome:R-HSA-5637704
#> 3180 reactome:R-HSA-5637694
#> 3181 reactome:R-HSA-5637693
#> 3182 reactome:R-HSA-5637698
#> 3183 reactome:R-HSA-5637706
#> 3184 reactome:R-HSA-1248675
#> 3185 reactome:R-HSA-1225979
#> 3395 reactome:R-HSA-9646489
#> 3396 reactome:R-HSA-9654429
#> 3397 reactome:R-HSA-9654438
#> 3399 reactome:R-HSA-9664564
#> 3400 reactome:R-HSA-9664561
#> 3401 reactome:R-HSA-9664562
#> 3403 reactome:R-HSA-9664609
#> 3404 reactome:R-HSA-9664598
#> 3405 reactome:R-HSA-9664604
#> 3406 reactome:R-HSA-9664613
#> 3407 reactome:R-HSA-9664603
#> 3408 reactome:R-HSA-9664614
#> 3417 reactome:R-HSA-9664917
#> 3418 reactome:R-HSA-9664921
#> 3419 reactome:R-HSA-9664932
#> 3423 reactome:R-HSA-9665001
#> 3424 reactome:R-HSA-9664996
#> 3425 reactome:R-HSA-9665037
#> 3426 reactome:R-HSA-9665033
#> 3438 reactome:R-HSA-9665391
#> 3439 reactome:R-HSA-9665380
#> 3442 reactome:R-HSA-9665436
#> 3443 reactome:R-HSA-9665441
#> 3444 reactome:R-HSA-9665437
#> 3445 reactome:R-HSA-9665449
#> 3446 reactome:R-HSA-9665446
#> 3447 reactome:R-HSA-9665453
#> 3448 reactome:R-HSA-9665460
#> 3450 reactome:R-HSA-9665867
#> 3451 reactome:R-HSA-9665877
#> 3452 reactome:R-HSA-9665910
#> 3453 reactome:R-HSA-9665887
#> 3454 reactome:R-HSA-9665882
#> 3455 reactome:R-HSA-9665869
#> 3456 reactome:R-HSA-9665917
#> 3457 reactome:R-HSA-9665916
#> 3458 reactome:R-HSA-9665925
#> 3459 reactome:R-HSA-9665924
#>                                                                                                                                                name
#> 2                                                                                                                                         EGFR:AAMP
#> 3                                                                                                                             EGFR:EGF-like ligands
#> 4                                                                                                                       EGF-like ligands:EGFR dimer
#> 5                                                                                                                  EGF-like ligands:p-6Y EGFR dimer
#> 6                                                                                                                         Activated EGFR:PLC-gamma1
#> 7                                                                                                                 Activated EGFR:Phospho-PLC-gamma1
#> 9                                                                                                                   GRB2:SOS-EGF-Phospho-EGFR dimer
#> 12                                                                                                                  EGF-like ligands:p-6Y EGFR:SHC1
#> 13                                                                                                                     EGF:Phospho-EGFR:Phospho-SHC
#> 14                                                                                                      GRB2:SOS-Phospho-SHC:EGF:Phospho-EGFR dimer
#> 18                                                                                                                 GAB1:GRB2-EGF-Phospho-EGFR dimer
#> 19                                                                                                         GRB2:Phospho GAB1-EGF-Phospho-EGFR dimer
#> 24                                                                                                   EGF-like ligands:p-6Y-EGFR:GRB2:p-5Y-GAB1:PI3K
#> 25                                                                                                    SHP2-GRB2:Phospho GAB1-EGF-Phospho-EGFR dimer
#> 26                                                                                             EGF-like ligands:p-6Y-EGFR:GRB2:p-Y627,659-GAB1:SHP2
#> 30                                                                                                   EGF-like ligands:p-5Y-EGFR:GRB2:p-5Y-GAB1:SHP2
#> 31                                                                                                      EGF-like ligands:Phospho-EGFR (-Y992) dimer
#> 32                                                                                                           EGF-like ligands:Phospho-EGFR (- Y992)
#> 35                                                                                                                   EGF-like ligands:p-6Y-EGFR:CBL
#> 36                                                                                                            EGF-like ligands:p-6Y-EGFR:p-Y371-CBL
#> 37                                                                                               EGF-like ligands:p-6Y-EGFR:CBL:Ub-p-Y53/55-SPRY1/2
#> 39                                                                                                         EGF-like ligands:Ub-p-6Y-EGFR:p-Y371-CBL
#> 41                                                                                                              EGF-like ligands:p-6Y-EGFR:CBL:GRB2
#> 42                                                                                                       EGF-like ligands:p-6Y-EGFR:p-Y371-CBL:GRB2
#> 43                                                                                                    EGF-like ligands:Ub-p-6Y-EGFR:p-Y371-CBL:GRB2
#> 49                                                                                      EGF-like ligands:p-6Y-EGFR:p-Y371-CBL:GRB2:CIN85:Endophilin
#> 56                                                                                                       EGF:p-6Y-EGFR:CBL:Beta-Pix:CDC42:GTP:CIN85
#> 57                                                                                                             EGF:p-6Y-EGFR:CBL:Beta-Pix:CDC42:GTP
#> 60                                                                       EGF-like ligands:p-6Y-EGFR:p-Y371-CBL:CIN85:Endophilin:Epsin:Eps15L1:Eps15
#> 61                                                                                                       EGF:Phospho-EGFR dimer:CBL:Phospho-Sprouty
#> 62                                                                    EGF-like ligands:p-6Y-EGFR:p-Y371-CBL:Ub-CIN85:Endophilin:Epsin:Eps15L1:Eps15
#> 65                                                               EGF-like ligands:p-6Y-EGFR:p-Y371-CBL:CIN85:SPRY1/2:Endophilin:Epsin:Eps15L1:EPS15
#> 66                                                                                                                              EGF:p-6Y-EGFR dimer
#> 67                                                                                                                                 EGF:Phospho-EGFR
#> 68                                                                                               EGF:p-EGFR dimer dephosphorylated at Y1172 (Y1148)
#> 69                                                                                                          EGF:p-Y992,Y1045,Y1068,Y1086,Y1173-EGFR
#> 70                                                          EGF-like ligands:p-6Y-EGFR:p-Y371-CBL:GRB2:CIN85:Endophilin:EPN1:EPS15L1:EPS15:HGS:STAM
#> 71                                                        EGF-like ligands:p-6Y-EGFR:p-Y371-CBL:GRB2:CIN85:Endophilin:EPN1:EPS15L1:p-EPS15:HGS:STAM
#> 531                                                                                                                                    VEGFR1 dimer
#> 534                                                                                                                                    VEGFR2 dimer
#> 537                                                                                                                                    VEGFR3 dimer
#> 539                                                                                                                          VEGFA-165 dimer:VEGFR2
#> 541                                                                                                                   VEGFA dimer:p-6Y-VEGFR2 dimer
#> 542                                                                                                                               p-6Y-VEGFR2 dimer
#> 543                                                                                                               VEGFA dimer:p-6Y-VEGFR2 dimer:SHB
#> 544                                                                                                           VEGFA dimer:p-6Y-VEGFR2 dimer:p-S-SHB
#> 545                                                                                                                    VEGFA:p-6Y-VEGFR2:p-SHB:PTK2
#> 546                                                                                                             VEGFA:p-6Y-VEGFR2:p-SHB:p-Y397-PTK2
#> 547                                                                                                       VEGFA:p-6Y-VEGFR2:p-SHB:p-Y397-PTK2:SRC-1
#> 548                                                                                                         VEGFA:p-6Y-VEGFR2:p-SHB:p-5Y-PTK2:SRC-1
#> 549                                                                                                               VEGFA:p-6Y-VEGFR2:p-SHB:p-5Y-PTK2
#> 551                                                                                                          VEGFA:p-6Y-VEGFR2:Integrin alphaVbeta3
#> 552                                                                                                VEGFA:p-6Y-VEGFR2:p-SHB:p-5Y-PTK2:SRC-1:HSP90AA1
#> 557                                                                                          VEGFA:p-6Y-VEGFR2:pS-SHB:p-5Y,S732-PTK2:SRC-1:HSP90AA1
#> 558                                                                                                         VEGFA:p-6Y-VEGFR2:pS-SHB:p-5Y,S732-PTK2
#> 559                                                                                                    VEGFA:p-6Y-VEGFR2:Integrin alphaVbeta3:PTK2B
#> 560                                                                                             VEGFA:p-6Y-VEGFR2:Integrin alphaVbeta3:p-Y402-PTK2B
#> 561                                                                                           VEGFA:p-6Y-VEGFR2:p-SHB:p-6Y,S732-PTK2:SRC-1:HSP90AA1
#> 562                                                                                                          VEGFA:p-6Y-VEGFR2:p-SHB:p-6Y,S732-PTK2
#> 563                                                                                               VEGFA:p-6Y-VEGFR2:p-SHB:p-7Y-PTK2:SRC-1:HSP90:PXN
#> 564                                                                                    VEGFA:p-6Y-VEGFR2:p-SHB:p-7Y-PTK2:SRC-1:HSP90:p-Y31,Y118-PXN
#> 565                                                                                     VEGFA:p-6Y-VEGFR2:p-SHB:p-6Y,S732-PTK2:SRC-1:HSP90AA1:BCAR1
#> 566                                                                                       VEGFA:p-6Y-VEGFR2:p-SHB:p-7Y-PTK2:SRC-1:HSP90:p-12Y-BCAR1
#> 567                          VEGFA:p-6Y-VEGFR2:p-SHB:p-7Y-PTK2:SRC-1:HSP90:p-12Y-BCAR1,VEGFA:p-6Y-VEGFR2:p-SHB:p-7Y-PTK2:SRC-1:HSP90:p-Y31,Y118-PXN
#> 568                                                                    VEGFA:p-6Y-VEGFR2:p-SHB:p-7Y-PTK2:SRC-1:HSP90:p-12Y-BCAR1/p-Y31,Y118-PAX:CRK
#> 570  CRK:DOCK180:ELMO1,ELMO2:VEGFA:p-6Y-VEGFR2:p-SHB:p-7Y-PTK2:SRC-1:HSP90:p-12Y-BCAR1,VEGFA:p-6Y-VEGFR2:p-SHB:p-7Y-PTK2:SRC-1:HSP90:p-Y31,Y118-PXN
#> 574                                                                                                                          VEGFA:p-6Y-VEGFR2:SHC2
#> 575                                                                                                            VEGFA dimer:p-6Y-VEGFR2 dimer:SH2D2A
#> 576                                                                                                      VEGFA dimer:p-6Y-VEGFR2 dimer:SH2D2A:SRC-1
#> 577                                                                                                           VEGFA:p-6Y-VEGFR2:SH2D2A:p-Y418-SRC-1
#> 578                                                                                                       VEGFA:p-6Y-VEGFR2:SH2D2A:p-Y418-SRC-1:AXL
#> 579                                                                                           VEGFA:p-6Y-VEGFR2:SH2D2A:p-Y418-SRC-1:p-Y779,Y821-AXL
#> 580                                                                                      VEGFA:p-6Y-VEGFR2:SH2D2A:p-Y418-SRC-1:p-Y779,Y821-AXL:PI3K
#> 581                                                                                                              VEGFA dimer:p-6Y-VEGFR2 dimer:PI3K
#> 615                                                                                                         VEGFA dimer:p-6Y-VEGFR2 dimer:NCK1,NCK2
#> 616                                                                                                     VEGFA dimer:p-6Y-VEGFR2 dimer:NCK1,NCK2:FYN
#> 617                                                                                              VEGFA dimer:p-6Y-VEGFR2 dimer:NCK1,NCK2:p-Y420-FYN
#> 618                                                                                          VEGFA dimer:p-6Y-VEGFR2 dimer:NCK1,NCK2:p-S21,Y420-FYN
#> 619                                                                                     VEGFA dimer:p-6Y-VEGFR2 dimer:NCK1,NCK2:p-S21,Y420-FYN:PAK2
#> 622                                                                                             VEGFA:p-6Y-VEGFR2:NCK:p-S21,Y420-FYN:PAK2:CDC42:GTP
#> 624                                                                                        VEGFA:p-6Y-VEGFR2:NCK:p-S21,Y420-FYN:p-2Y-PAK2:CDC42:GTP
#> 627                                                                                        VEGFA:p-6Y-VEGFR2:NCK:p-S21,Y420-FYN:p-3Y-PAK2:CDC42:GTP
#> 631                                                                                                             VEGFA dimer:p-6Y-VEGFR2 dimer:PLCG1
#> 632                                                                                                        VEGFA dimer:p-6Y-VEGFR2 dimer:p-4Y-PLCG1
#> 643                                                                                                                               NRP1:VEGFR2 dimer
#> 644                                                                                                                               NRP2:VEGFR1 dimer
#> 688                                                                                                               Ligand-Activated EGFR/ERBB3/ERBB4
#> 689                                                                                                                                        EGF:EGFR
#> 695                                                                                                                                  EGF:EGFR:ERBB2
#> 701                                                                                                                           EGF:EGFR:p-Y877-ERBB2
#> 705                                                                                                                  EGF:p-6Y-EGFR:p-7Y,Y1112-ERBB2
#> 706                                                                                                                                   EGF:p-6Y-EGFR
#> 713                                                                                                                  EGF:p-6Y-EGFR:p-6Y,Y1112-ERBB2
#> 728                                                                                                          Phosphorylated ERBB2:EGFR heterodimers
#> 729                                                                                                                        GRB2:SOS1:P-ERBB2:P-EGFR
#> 736                                                                                                                    EGF:p-EGFR:p-ERBB2:GRB2:GAB1
#> 737                                                                                                                   PI3K:GRB2:GAB1:P-ERBB2:P-EGFR
#> 738                                                                                                                PI3Kp85:GRB2:GAB1:P-ERBB2:P-EGFR
#> 739                                                                                                                        EGF:p-EGFR:p-ERBB2:PLCG1
#> 749                                                                   EGF:p-Y992,Y1045,Y1068,Y1086,Y1173-EGFR:p-Y1023,Y1139,Y1196,Y1221,Y1222-ERBB2
#> 752                                                              EGF:p-Y992,Y1045,Y1068,Y1086,Y1173-EGFR:p-Y877,Y1023,Y1139,Y1196,Y1221,Y1222-ERBB2
#> 755                                                                                                   Phosphorylated ERBB2:EGFR heterodimers:PTPN18
#> 756                                                        Phosphorylated ERBB2(dephosphorylated at Y1196,Y1112 and Y1248):EGFR heterodimers:PTPN18
#> 757                                                               Phosphorylated ERBB2(dephosphorylated at Y1196,Y1112 and Y1248):EGFR heterodimers
#> 758                                                                                                   EGF:p-6Y-EGFR:p-Y1023,Y1139,Y1221,Y1222-ERBB2
#> 759                                                                                              EGF:p-6Y-EGFR:p-Y877,Y1023,Y1139,Y1221,Y1222-ERBB2
#> 778                                                                                                                          ERBB4:EGFR heterodimer
#> 1532                                                                                                                         HB-EGF:p-6Y-EGFR dimer
#> 1533                                                                                                                               HB-EGF:p-6Y-EGFR
#> 1534                                                                                                                     GRB2:SOS1:HB-EGF:p-6Y-EGFR
#> 1729                                                                                                                      EGF:p-6Y-EGFR dimer:NICD3
#> 2256                                                                                                                EGF-like ligands:p-6Y EGFR:PTK2
#> 2257                                                                                                   EGF-like ligands:p-6Y EGFR dimer:p-Y397 PTK2
#> 2941                                                                                                                                     HBEGF:EGFR
#> 2942                                                                                                                               HBEGF:EGFR:GPNMB
#> 2943                                                                                                                        HBEGF:EGFR:p-Y525-GPNMB
#> 2944                                                                                                   HBEGF:EGFR:p-Y525-GPNMB:LINC01139:PTK6:LRRK2
#> 2945                                                                                            HBEGF:EGFR:p-Y525-GPNMB:LINC01139:p-Y351-PTK6:LRRK2
#> 3135                                                                                                     Ligand-responsive EGFR mutants:HSP90:CDC37
#> 3136                                                                      Ligand-responsive EGFR mutants resistant to non-covalent TKIs:HSP90:CDC37
#> 3137                                                                      Ligand-responsive EGFR mutants sensitive to non-covalent TKIs:HSP90:CDC37
#> 3138                                                                                                 EGF:Ligand-responsive EGFR mutants:HSP90:CDC37
#> 3139                                                                  EGF:Ligand-responsive EGFR mutants resistant to non-covalent TKIs:HSP90:CDC37
#> 3140                                                                  EGF:Ligand-responsive EGFR mutants sensitive to non-covalent TKIs:HSP90:CDC37
#> 3141                                                                                                       EGF:Ligand-responsive EGFR mutants dimer
#> 3142                                                                                                           Ligand-responsive EGFR mutants dimer
#> 3143                                                                                                Active dimers of ligand-responsive EGFR mutants
#> 3144                                                                                                      Ligand-responsive p-6Y-EGFR mutant dimers
#> 3145                                                                                              Dimers of EGF:Ligand-responsive p-6Y-EGFR mutants
#> 3146                                                                                            EGF:Ligand-responsive p-6Y-EGFR mutants:HSP90:CDC37
#> 3147                                                             EGF:Ligand-responsive p-6Y-EGFR mutants resistant to non-covalent TKIs:HSP90:CDC37
#> 3148                                                             EGF:Ligand-responsive p-6Y-EGFR mutants sensitive to non-covalent TKIs:HSP90:CDC37
#> 3149                                                                                                  Dimers of ligand-responsive p-6Y-EGFR mutants
#> 3150                                                                                                Ligand-responsive p-6Y-EGFR mutants:HSP90:CDC37
#> 3151                                                                 Ligand-responsive p-6Y-EGFR mutants resistant to non-covalent TKIs:HSP90:CDC37
#> 3152                                                                 Ligand-responsive p-6Y-EGFR mutants sensitive to non-covalent TKIs:HSP90:CDC37
#> 3153                                                                                                        GRB2:SOS1:Phospho-EGFRmutKD:HSP90:CDC37
#> 3154                                                                                                       Ligand-responsive p-6Y-EGFR mutants:SHC1
#> 3155                                                                                            Ligand-responsive p-6Y-EGFR mutants:p-Y349,350-SHC1
#> 3156                                                                                           GRB2:SOS1:Phospho-SHC1:Phospho-EGFRmutKD:HSP90:CDC37
#> 3157                                                                                                  Ligand-responsive p-6Y-EGFR mutants:GRB2:GAB1
#> 3158                                                                                             Ligand-responsive p-6Y-EGFR mutants:GRB2:GAB1:PI3K
#> 3159                                                                                                                    Phospho-EGFRmutKD:PLCgamma1
#> 3160                                                                                  Ligand-responsive p-6Y-EGFR mutants:p-Y472,771,783,1254-PLCG1
#> 3161                                                                                                        Ligand-responsive p-6Y-EGFR mutants:CBL
#> 3162                                                                                                 Ligand-responsive p-6Y-EGFR mutants:p-Y371-CBL
#> 3164                                                                                                       Resistant ligand-responsive EGFR mutants
#> 3165                                                                    Dimers of EGF:Ligand-responsive EGFR mutants resistant to non-covalent TKIs
#> 3166                                                                        Dimers of ligand-responsive EGFR mutants resistant to non-covalent TKIs
#> 3167                                                                                Sensitive ligand-responsive EGFR mutants:Non-covalent EGFR TKIs
#> 3168                                                                 Active dimers of ligand-responsive EGFR mutants sensitive to non-covalent TKIs
#> 3169                                                                    Dimers of EGF:Ligand-responsive EGFR mutants sensitive to non-covalent TKIs
#> 3170                                                                        Dimers of ligand-responsive EGFR mutants sensitive to non-covalent TKIs
#> 3171                                                                                    Resistant ligand-responsive EGFR mutants:Covalent EGFR TKIs
#> 3172                                                                                                                        EGFRvIIImut:HSP90:CDC37
#> 3173                                                                                                                          EGFRvIII mutant dimer
#> 3174                                                                                                                     p-5Y-EGFRvIII mutant dimer
#> 3175                                                                                                                Phospho-EGFRvIIImut:HSP90:CDC37
#> 3176                                                                                                                        p-5Y-EGFRvIII:GRB2:SOS1
#> 3177                                                                                                                             p-5Y-EGFRvIII:SHC1
#> 3178                                                                                                                 p-5Y-EGFRvIII:p-Y349,Y350-SHC1
#> 3179                                                                                                       p-5Y-EGFRvIII:p-Y349,Y350-SHC1:GRB2:SOS1
#> 3180                                                                                                                        p-5Y-EGFRvIII:GRB2:GAB1
#> 3181                                                                                                                   p-5Y-EGFRvIII:GRB2:GAB1:PI3K
#> 3182                                                                                                                            p-5Y-EGFRvIII:PLCG1
#> 3183                                                                                                          p-5Y-EGFRvIII:p-Y771,Y783,Y1254-PLCG1
#> 3184                                                                                                                                 EGFR:Cetuximab
#> 3185                                                                                                                       EGF:EGFR Dimer:irrevTKIs
#> 3395                                                                                                            Ligand-Activated EGFR,ERBB3,(ERBB4)
#> 3396                                                                                                            Ligand-Activated EGFR,(ERBB3,ERBB4)
#> 3397                                                                                                            Ligand-Activated ERBB3,(EGFR,ERBB4)
#> 3399                                                                                           ERBB2 KD mutants:Ligand-Activated EGFR,ERBB3,(ERBB4)
#> 3400                                                                                           ERBB2 KD mutants:Ligand-Activated EGFR,(ERBB3,ERBB4)
#> 3401                                                                                           ERBB2 KD mutants:Ligand-Activated ERBB3,(EGFR,ERBB4)
#> 3403                                                                                p-6Y-ERBB2 KD mutants:Ligand-Activated p-EGFR,(p-ERBB3,p-ERBB4)
#> 3404                                                                                                      Ligand-Activated p-EGFR,(p-ERBB3,p-ERBB4)
#> 3405                                                                                p-6Y-ERBB2 KD mutants:Ligand-Activated p-EGFR,p-ERBB3,(p-ERBB4)
#> 3406                                                                                                      Ligand-Activated p-EGFR,p-ERBB3,(p-ERBB4)
#> 3407                                                                                p-6Y-ERBB2 KD mutants:Ligand-Activated p-ERBB3,(p-EGFR,p-ERBB4)
#> 3408                                                                                                      Ligand-Activated p-ERBB3,(p-EGFR,p-ERBB4)
#> 3417                                                                                                            p-6Y-ERBB2 KD mutants:EGF:p-6Y-EGFR
#> 3418                                                                                                  p-6Y-ERBB2 KD mutants:EGF:p-6Y-EGFR:GRB2:GAB1
#> 3419                                                                                             p-6Y-ERBB2 KD mutants:EGF:p-6Y-EGFR:GRB2:GAB1:PI3K
#> 3423                                                                                                            p-6Y-ERBB2 KD mutants:EGF:p-6Y-EGFR
#> 3424                                                                                                  p-6Y-ERBB2 KD mutants:EGF:p-6Y-EGFR:GRB2:SOS1
#> 3425                                                                                                    p-6Y-ERBB2 KD mutants (PLCG1):EGF:p-6Y-EGFR
#> 3426                                                                                                      p-6Y-ERBB2 KD mutants:EGF:p-6Y-EGFR:PLCG1
#> 3438                                                                                                                     ERBB2 ECD mutants:EGF:EGFR
#> 3439                                                                                                           p-6Y-ERBB2 ECD mutants:EGF:p-6Y-EGFR
#> 3442                                                                                                      p-6Y-ERBB2 ECD mutants:EGF:p-6Y-EGFR:SHC1
#> 3443                                                                                          p-6Y-ERBB2 ECD mutants:EGF:p-6Y-EGFR:p-Y349,Y350-SHC1
#> 3444                                                                                p-6Y-ERBB2 ECD mutants:EGF:p-6Y-EGFR:p-Y349,Y350-SHC1:GRB2:SOS1
#> 3445                                                                                                 p-6Y-ERBB2 ECD mutants:EGF:p-6Y-EGFR:GRB2:SOS1
#> 3446                                                                                                 p-6Y-ERBB2 ECD mutants:EGF:p-6Y-EGFR:GRB2:GAB1
#> 3447                                                                                            p-6Y-ERBB2 ECD mutants:EGF:p-6Y-EGFR:GRB2:GAB1:PI3K
#> 3448                                                                                                     p-6Y-ERBB2 ECD mutants:EGF:p-6Y-EGFR:PLCG1
#> 3450                                                                                      ERBB2 TMD/JMD mutants:Ligand-Activated EGFR,ERBB3,(ERBB4)
#> 3451                                                                           p-6Y-ERBB2 TMD/JMD mutants:Ligand-Activated p-EGFR,p-ERBB3,(p-ERBB4)
#> 3452                                                                     p-6Y-ERBB2 TMD/JMD mutants (RAS):Ligand-Activated p-EGFR,p-ERBB3,(p-ERBB4)
#> 3453                                                                p-6Y-ERBB2 TMD/JMD mutants (RAS):Ligand-Activated p-EGFR,p-ERBB3,(p-ERBB4):SHC1
#> 3454                                                                    p-6Y-ERBB2 TMD/JMD mutants:Ligand-Activated p-EGFR,p-ERBB3,(p-ERBB4):p-SHC1
#> 3455                                                          p-6Y-ERBB2 TMD/JMD mutants:Ligand-Activated p-EGFR,p-ERBB3,(p-ERBB4):p-SHC1:GRB2:SOS1
#> 3456                                                                                                 p-6Y-ERBB2 TMD/JMD mutants (RAS):EGF:p-6Y-EGFR
#> 3457                                                                                       p-6Y-ERBB2 TMD/JMD mutants (RAS):EGF:p-6Y-EGFR:GRB2:SOS1
#> 3458                                                                                                                 p-6Y-ERBB2 R678Q:EGF:p-6Y-EGFR
#> 3459                                                                                                           p-6Y-ERBB2 R678Q:EGF:p-6Y-EGFR:PLCG1
```
