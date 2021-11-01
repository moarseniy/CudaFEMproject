*KEYWORD
$
$ File created by cubit version 17.1
$...>....1....>....2....>....3....>....4....>....5....>....6....>....7....>....8
*TITLE
$ cubit(D:/Fidesys/result1.k): 03/14/2021: 18:44:06
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
      26               5               0               0       0       0
      21               0               5               0       0       0
      27               0        -1.66667               0       0       0
      23              -5               0               0       0       0
      20              -5               5               0       0       0
      25               0              -5               0       0       0
      22              -5              -5               0       0       0
      24               5              -5               0       0       0
      19               5               5               0       0       0
$
$...>....1....>....2....>....3....>....4....>....5....>....6....>....7....>....8
$    eid     pid      n1      n2      n3      n4      n5      n6      n7      n8
*ELEMENT_SHELL
       1       1      26      21      27      27
       2       1      21      23      27      27
       3       1      21      20      23      23
       4       1      25      27      22      22
       5       1      26      27      24      24
       6       1      26      19      21      21
       7       1      27      25      24      24
       8       1      27      23      22      22
*BOUNDARY_SPC_SET
         1         0         1         1         0         0         0         0
*SET_NODE_LIST
         1
        22
*BOUNDARY_SPC_SET
         2         0         0         1         0         0         0         0
*SET_NODE_LIST
         2
        24
*LOAD_NODE_SET
         3         2         1       1.0
*DEFINE_CURVE
         1
       0.0         1
*SET_NODE_LIST
         3
        20        19        21
*END
