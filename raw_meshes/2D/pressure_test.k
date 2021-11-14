*KEYWORD
$
$ File created by cubit version 17.1
$...>....1....>....2....>....3....>....4....>....5....>....6....>....7....>....8
*TITLE
$ cubit(D:/Fidesys/pressure_test/pressure_test.k): 11/14/2021: 23:04:58
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
      15              -5     8.88178e-16               0       0       0
      12              -5              -5               0       0       0
      16    -8.88178e-16              -5               0       0       0
      11              -5               5               0       0       0
      18        -0.46426         0.46426               0       0       0
      14     8.88178e-16               5               0       0       0
      17               5    -8.88178e-16               0       0       0
      10               5               5               0       0       0
      13               5              -5               0       0       0
$
$...>....1....>....2....>....3....>....4....>....5....>....6....>....7....>....8
$    eid     pid      n1      n2      n3      n4      n5      n6      n7      n8
*ELEMENT_SHELL
       1       1      15      12      16      16
       2       1      11      18      14      14
       3       1      17      10      14      14
       4       1      11      15      18      18
       5       1      13      17      16      16
       6       1      17      18      16      16
       7       1      18      15      16      16
       8       1      18      17      14      14
*LOAD_SEGMENT
         1       1.0       0.0        11        15        15        15       0.0
*DEFINE_CURVE
         1
       0.0     1e+09
*LOAD_SEGMENT
         2       1.0       0.0        15        12        12        12       0.0
*DEFINE_CURVE
         2
       0.0     1e+09
*LOAD_SEGMENT
         3       1.0       0.0        13        17        17        17       0.0
*DEFINE_CURVE
         3
       0.0     1e+09
*LOAD_SEGMENT
         4       1.0       0.0        17        10        10        10       0.0
*DEFINE_CURVE
         4
       0.0     1e+09
*BOUNDARY_SPC_SET
         1         0         1         1         0         0         0         0
*SET_NODE_LIST
         1
        10        11        12        13        14        16
*END
