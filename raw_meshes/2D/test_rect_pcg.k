*KEYWORD
$
$ File created by cubit version 2021.4
$...>....1....>....2....>....3....>....4....>....5....>....6....>....7....>....8
*TITLE
$ cubit(C:/Users/mokin/fidesys/test_rect_pcg/test_rect_pcg.k): 03/24/2022: 02:13:2
*PART
Block 1
$    pid     sid     mid   eosid    hgid    grav  adpopt    tmid
       1       1     100       0       0       0       0     100
*SECTION_SHELL
       1
$     t1      t2      t3      t4
       1       1       1       1

*MAT_ELASTIC
     100   7e-06  206800    0.29     0.0     0.0       0
*MAT_THERMAL_ISOTROPIC
     100     0.0     0.0     0.0     0.0     0.0
     500   0.045
*NODE
$    nid               x               y               z      tc      rc
     385              -5             2.5               0       0       0
     389    -1.57182e-07            -2.5               0       0       0
     388     1.57182e-07             2.5               0       0       0
     387               5            -2.5               0       0       0
     384               5             2.5               0       0       0
     386              -5            -2.5               0       0       0
$
$...>....1....>....2....>....3....>....4....>....5....>....6....>....7....>....8
$    eid     pid      n1      n2      n3      n4      n5      n6      n7      n8
*ELEMENT_SHELL
       1       1     385     389     388     388
       2       1     387     384     388     388
       3       1     389     385     386     386
       4       1     387     388     389     389
*LOAD_SEGMENT
         1       1.0       0.0       384       388       388       388       0.0
*DEFINE_CURVE
         1
       0.0      7000
*LOAD_SEGMENT
         2       1.0       0.0       384       388       388       388       0.0
*DEFINE_CURVE
         2
       0.0      7000
*BOUNDARY_SPC_SET
         1         0         1         1         0         0         0         0
*SET_NODE_LIST
         1
       385       386
*END
