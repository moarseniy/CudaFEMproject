*KEYWORD
$
$ File created by cubit version 17.1
$...>....1....>....2....>....3....>....4....>....5....>....6....>....7....>....8
*TITLE
$ cubit(D:/Fidesys/small_rect/small_rect.k): 03/31/2021: 22:40:49
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
     144              -5     8.88178e-16               0       0       0
     141              -5              -5               0       0       0
     145    -8.88178e-16              -5               0       0       0
     140              -5               5               0       0       0
     147        -0.46426         0.46426               0       0       0
     143     8.88178e-16               5               0       0       0
     146               5    -8.88178e-16               0       0       0
     139               5               5               0       0       0
     142               5              -5               0       0       0
$
$...>....1....>....2....>....3....>....4....>....5....>....6....>....7....>....8
$    eid     pid      n1      n2      n3      n4      n5      n6      n7      n8
*ELEMENT_SHELL
       1       1     144     141     145     145
       2       1     140     147     143     143
       3       1     146     139     143     143
       4       1     140     144     147     147
       5       1     142     146     145     145
       6       1     146     147     145     145
       7       1     147     144     145     145
       8       1     147     146     143     143
*BOUNDARY_SPC_SET
         1         0         1         1         0         0         0         0
*SET_NODE_LIST
         1
       139       143       140
*LOAD_NODE_SET
         2         2         1       1.0
*DEFINE_CURVE
         1
       0.0        -1
*SET_NODE_LIST
         2
       141       145       142
*END
