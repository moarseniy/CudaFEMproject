#include "tests.h"

// 18 global DOFs, 12 local DOFs
// pass void elems, m_e, b_e
void test_EbePCG_diag(std::vector<Element> &elems, std::vector<float*> &m_e, std::vector<float*> &b_e, float soln[18])
{
  Element elem1, elem2, elem3;

  elem1.nodesIds[0] = 6;
  elem1.nodesIds[1] = 9;
  elem1.nodesIds[2] = 0;
  elem1.nodesIds[3] = 3;

  elem2.nodesIds[0] = 15;
  elem2.nodesIds[1] = 12;
  elem2.nodesIds[2] = 6;
  elem2.nodesIds[3] = 3;

  elem3.nodesIds[0] = 6;
  elem3.nodesIds[1] = 12;
  elem3.nodesIds[2] = 3;
  elem3.nodesIds[3] = 9;

  // Specify the diagonal elements
  float diag1[12] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0 };
  float diag2[12] = { 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0 };
  float diag3[12] = { 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0 };
  for (int i = 0; i < 12; ++i) elem1.Klocal.Set(i,i,diag1[i]);
  for (int i = 0; i < 12; ++i) elem2.Klocal.Set(i,i,diag2[i]);
  for (int i = 0; i < 12; ++i) elem3.Klocal.Set(i,i,diag3[i]);
  elems.push_back(elem1); elems.push_back(elem2); elems.push_back(elem3);
  // Specify right-hand vector b_e
  //int *array{ new int[5]{ 10, 7, 15, 3, 11 } };
  float *b_1{ new float[12]{ 4.24, 8.2, -0.13, 0.0, 5.03, -0.04, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 }};
  float *b_2{ new float[12]{ 12.3, 0.0, 0.578, 3.0, 3.01, 0.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0 }};
  float *b_3{ new float[12]{ 0.0, 0.0, -2.89, 0.0, 9.15, 0.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0 }};
  std::cout << b_1[0] << std::endl;
  b_e.push_back(b_1); b_e.push_back(b_2); b_e.push_back(b_3);

  float m[18];    // Diagonal elements of the global matrix. In this case, it's the global matrix itself.

  m[0]    = elem1.Klocal.Get(6,6);
  m[1]    = elem1.Klocal.Get(7,7);
  m[2]    = elem1.Klocal.Get(8,8);

  m[3]    = elem1.Klocal.Get(9,9)    + elem2.Klocal.Get(9,9)     + elem3.Klocal.Get(6,6);
  m[4]    = elem1.Klocal.Get(10,10)  + elem2.Klocal.Get(10,10)   + elem3.Klocal.Get(7,7);
  m[5]    = elem1.Klocal.Get(11,11)  + elem2.Klocal.Get(11,11)   + elem3.Klocal.Get(8,8);

  m[6]    = elem1.Klocal.Get(0,0)    + elem2.Klocal.Get(6,6)     + elem3.Klocal.Get(0,0);
  m[7]    = elem1.Klocal.Get(1,1)    + elem2.Klocal.Get(7,7)     + elem3.Klocal.Get(1,1);
  m[8]    = elem1.Klocal.Get(2,2)    + elem2.Klocal.Get(8,8)     + elem3.Klocal.Get(2,2);

  m[9]    = elem1.Klocal.Get(3,3)                                + elem3.Klocal.Get(9,9);
  m[10]   = elem1.Klocal.Get(4,4)                                + elem3.Klocal.Get(10,10);
  m[11]   = elem1.Klocal.Get(5,5)                                + elem3.Klocal.Get(11,11);

  m[12]   =                            elem2.Klocal.Get(3,3)     + elem3.Klocal.Get(3,3);
  m[13]   =                            elem2.Klocal.Get(4,4)     + elem3.Klocal.Get(4,4);
  m[14]   =                            elem2.Klocal.Get(5,5)     + elem3.Klocal.Get(5,5);

  m[15]   =                            elem2.Klocal.Get(0,0);
  m[16]   =                            elem2.Klocal.Get(1,1);
  m[17]   =                            elem2.Klocal.Get(2,2);

  // Think how to get it without finding global matrix m[18].
  //float m_1[12], m_2[12], m_3[12];
  float *m_1 = new float [12];
  float *m_2 = new float [12];
  float *m_3 = new float [12];
  for (int i = 0; i < 4; ++i) {
    for (int idim = 0; idim < 3; ++idim) {
      m_1[3*i + idim] = m[elem1.nodesIds[i] + idim];
      //std::cout << "m_1[" << 3*i + idim << "] = " << m_1[3*i + idim] << std::endl;
      m_2[3*i + idim] = m[elem2.nodesIds[i] + idim];
      //std::cout << "m_2[" << 3*i + idim << "] = " << m_2[3*i + idim] << std::endl;
      m_3[3*i + idim] = m[elem3.nodesIds[i] + idim];
      //std::cout << "m_3[" << 3*i + idim << "] = " << m_3[3*i + idim] << std::endl;
    }
  }
  m_e.push_back(m_1); m_e.push_back(m_2); m_e.push_back(m_3);

  float b[18];    // Right-hand vector of the global system.

  b[0]    = b_1[6];
  b[1]    = b_1[7];
  b[2]    = b_1[8];

  b[3]    = b_1[9]   + b_2[9]    + b_3[6];
  b[4]    = b_1[10]  + b_2[10]   + b_3[7];
  b[5]    = b_1[11]  + b_2[11]   + b_3[8];

  b[6]    = b_1[0]   + b_2[6]     + b_3[0];
  b[7]    = b_1[1]   + b_2[7]     + b_3[1];
  b[8]    = b_1[2]   + b_2[8]     + b_3[2];

  b[9]    = b_1[3]                + b_3[9];
  b[10]   = b_1[4]                + b_3[10];
  b[11]   = b_1[5]                + b_3[11];

  b[12]   =            b_2[3]     + b_3[3];
  b[13]   =            b_2[4]     + b_3[4];
  b[14]   =            b_2[5]     + b_3[5];

  b[15]   =            b_2[0];
  b[16]   =            b_2[1];
  b[17]   =            b_2[2];

  for (int i = 0; i < 18; ++i)
    soln[i] = m[i]/b[i];

  //delete [] m_1;
  //delete [] m_2;
  //delete [] m_3;
  //delete [] b_1;
  //delete [] b_2;
  //delete [] b_3;

}
