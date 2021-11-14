*KEYWORD
$
$ File created by cubit version 17.1
$...>....1....>....2....>....3....>....4....>....5....>....6....>....7....>....8
*TITLE
$ cubit(D:/Fidesys/udali_menya/udali_menya.k): 11/14/2021: 23:49:46
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
*LOAD_SEGMENT
         1       1.0       0.0       140       144       144       144       0.0
*DEFINE_CURVE
         1
       0.0     1e+09
*LOAD_SEGMENT
         2       1.0       0.0       144       141       141       141       0.0
*DEFINE_CURVE
         2
       0.0     1e+09
*LOAD_SEGMENT
         3       1.0       0.0       142       146       146       146       0.0
*DEFINE_CURVE
         3
       0.0     1e+09
*LOAD_SEGMENT
         4       1.0       0.0       146       139       139       139       0.0
*DEFINE_CURVE
         4
       0.0     1e+09
*BOUNDARY_SPC_SET
         1         0         1         1         0         0         0         0
*SET_NODE_LIST
         1
       143       145
*END
