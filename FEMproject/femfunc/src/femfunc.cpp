
#include "femfunc.h"

using namespace std;


float SetConstraints(int i, int j, float v, int index) {
  if (i == index || j == index) {
    return i == j ? 1.0f : 0.0f;
  } else {
    return v;
  }
}

void ApplyLoads(MyArray& F, const std::vector<Load>& loads) {
  CheckRunTime(__func__)
  for (std::vector<Load>::const_iterator it = loads.begin(); it != loads.end(); ++it) {
    F[it->dof] += it->value;
  }
}

void AssignLoadElement(FEMdataKeeper &FEMdata, std::unordered_map <int, std::vector<int>> nodeAdjElem) {
  CheckRunTime(__func__)
  for (std::vector<Load>::iterator it = FEMdata.loads.begin(); it != FEMdata.loads.end(); ++it) {
    it->assignElement(FEMdata.DIM, nodeAdjElem);
  }
}

// Make sure that elements are assigned to each load before calling this function!
// (call AssignLoadElement before calling this function)
// And make sure that loadVectors is empty before calling this function!
void GetMapElement2Loadvector(FEMdataKeeper &FEMdata, std::unordered_map <int, MyArray> &loadVectors, float t) {
  CheckRunTime(__func__)
  int DIM = FEMdata.DIM;
  for (std::vector<Load>::iterator it = FEMdata.loads.begin(); it != FEMdata.loads.end(); ++it) {
    assert(it->elem != -1);
    it->TimeDependentEntity::update(t);
    MyArray elem_load(6 * (DIM - 1));
    elem_load[ FEMdata.elements[it->elem].Global2LocalNode(it->dof / DIM) * DIM + it->dof % DIM ] = it->value;
    if (loadVectors.find(it->elem) == loadVectors.end())        // Key not present
      loadVectors[it->elem] = elem_load;
    else
      loadVectors[it->elem].add(elem_load);
  }
}

// Make sure that elements are assigned to each load before calling this function!
// (call AssignLoadElement before calling this function)
void ApplyLoads_EbE(FEMdataKeeper &FEMdata) {
  CheckRunTime(__func__)
  int DIM = FEMdata.DIM;
  for (std::vector<Load>::iterator it = FEMdata.loads.begin(); it != FEMdata.loads.end(); ++it) {
    //        int n_elems = nodeAdjElem[it->node].size();
    //        for (std::vector<int>::const_iterator it_elem = nodeAdjElem[it->node].begin(); it_elem != nodeAdjElem[it->node].end(); ++it_elem) {
    //            Element *elem = &FEMdata.elements[*it_elem];
    //            std::cout << *it_elem << std::endl;
    //            elem->Flocal[elem->Global2LocalNode(it->node) * 2 + 0] += it->x / static_cast<float>(n_elems);
    //            elem->Flocal[elem->Global2LocalNode(it->node) * 2 + 1] += it->y / static_cast<float>(n_elems);
    //        }
    assert(it->elem != -1);
    FEMdata.elements[it->elem].blocal[ FEMdata.elements[it->elem].Global2LocalNode(it->dof / DIM) * DIM + it->dof % DIM ] += it->value;
  }
}

void ApplyConstraints(SparseMatrixCOO& K, MyArray& F, const std::vector<Constraint>& constraints) {
  CheckRunTime(__func__)
  std::vector<int> indicesToConstraint;

  for (std::vector<Constraint>::const_iterator it = constraints.begin(); it != constraints.end(); ++it) {
    if (it->type & Constraint::UX) {
      indicesToConstraint.push_back(2 * it->node + 0);
    }
    if (it->type & Constraint::UY) {
      indicesToConstraint.push_back(2 * it->node + 1);
    }
  }


  //    for (int i = 0; i < indicesToConstraint.size(); i++) {
  //        cout << indicesToConstraint[i] << " ";
  //    }


  //    for (int i = 0; i < K.get_size(); i++) {
  //        for (int j = 0; j < indicesToConstraint.size(); j++) {
  //            if (K.get_x(i) == indicesToConstraint[j] || K.get_y(i) == indicesToConstraint[j]) {
  //                if (K.get_x(i) == K.get_y(i)) {
  //                    K.set_value(K.get_x(i), K.get_y(i), 1.0);
  //                } else {
  //                    K.set_value(K.get_x(i), K.get_y(i), 0.0);
  //                }
  //                //cout << K.get_x(i) << " " << K.get_y(i) << "\n\n";
  //            }
  //        }
  //    }

  int k = 0;
  int i = 0;
  int threshhold = K.get_size();
  while (i < threshhold) {
    //for (int i = 0; i < K.get_size(); i++) {
    for (int j = 0; j < indicesToConstraint.size(); j++) {
      if (K.get_x(i) == indicesToConstraint[j] || K.get_y(i) == indicesToConstraint[j]) {
        if (K.get_x(i) == K.get_y(i)) {
          //K.set_value(K.get_x(i), K.get_y(i), 1.0);
          k++;
        } else {
          //K.set_value(K.get_x(i), K.get_y(i), 0.0);
          K.pop(K.get_x(i), K.get_y(i));
          --threshhold; --i;
        }
        //cout << K.get_x(i) << " " << K.get_y(i) << "\n\n";
      }
    }

    ++i;
  }

  //for (int i = 0; i < n; ++i) {
  for (int j = 0; j < indicesToConstraint.size(); ++j) {
    //if (i == indicesToConstraint[j]) {
    F[indicesToConstraint[j]] = 0.0;
    //F[2 * indicesToConstraint[j] + 1] = 0.0;
    //std::cout << indicesToConstraint[j] << " ";
    //}
  }
  //std::cout << "num=" << indicesToConstraint.size() << "\n";
  //}
  //cout << "!!!!" << k << " " << indicesToConstraint.size() << " ";
}

SparseMatrixCOO AssemblyStiffnessMatrix(FEMdataKeeper &FEMdata) {
  CheckRunTime(__func__)
  vecSparseMatrixCOO triplets;
  for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        for (int ilocal = 0; ilocal < 2; ++ilocal) {
          for (int jlocal = 0; jlocal < 2; ++jlocal) {
            float value = it->Klocal(2 * i + ilocal, 2 * j + jlocal);
            if (value != 0.0) {
              Triplet tmp(2 * it->nodesIds[i] + ilocal, 2 * it->nodesIds[j] + jlocal, value);
              triplets.push_back(tmp);
            }
          }
        }
      }
    }
  }

  cout << "CalculateStiffnessMatrix success\n";
  cout << "Triplets Size = " << triplets.size() << std::endl;

  SparseMatrixCOO globalK(triplets.size());
  globalK.ConvertTripletToSparse(triplets);

  std::cout << "new size= "<<globalK.get_size()<<"\n";

  return globalK;
}

MyArray AssemblyF(FEMdataKeeper &FEMdata) {
  MyArray F(FEMdata.DIM * FEMdata.nodesCount);
  for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
    for (int i = 0; i < 3; ++i) {
      F[2 * it->nodesIds[i] + 0] += it->Flocal[2 * i + 0];
      F[2 * it->nodesIds[i] + 1] += it->Flocal[2 * i + 1];
    }
  }
  return F;
}

void AssemblyX(FEMdataKeeper &FEMdata, std::unordered_map <int, std::vector<int>> nodeAdjElem) {
  for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
    for (int i = 0; i < 3; ++i) {
      FEMdata.displacements[2 * it->nodesIds[i] + 0] += it->x[2 * i + 0];
      FEMdata.displacements[2 * it->nodesIds[i] + 1] += it->x[2 * i + 1];
    }
  }

  for (int global_node_num = 0; global_node_num < FEMdata.nodesCount; ++global_node_num) {
    float num_of_elems = static_cast <float> (nodeAdjElem[global_node_num].size());
    FEMdata.displacements[2 * global_node_num + 0] /= num_of_elems;
    FEMdata.displacements[2 * global_node_num + 1] /= num_of_elems;
  }
}

void AssemblyB(FEMdataKeeper &FEMdata, MyArray &b) {
  for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
    for (int i = 0; i < 3; ++i) {
      b[2 * it->nodesIds[i] + 0] += it->b[2 * i + 0];
      b[2 * it->nodesIds[i] + 1] += it->b[2 * i + 1];
    }
  }
}

void CalculateStressAndDeformation(int DIM,
                                   std::vector<MyArray> &Deformation,
                                   std::vector<MyArray> &Stress,
                                   std::vector<float> &epsilon_mises,
                                   std::vector<float> &sigma_mises,
                                   Matrix &D,
                                   std::vector<Element> &elements,
                                   MyArray &displacements, MyArray &all_B) {
  CheckRunTime(__func__)

  MyArray StressVector(3 * (DIM - 1));
  MyArray DeformationVector(3 * (DIM - 1));
  MyArray delta(6 * (DIM - 1));

  Matrix B(3 * (DIM - 1), 6 * (DIM - 1));

  for (int i = 0; i < elements.size(); ++i) {
    if (DIM == 2) {
      delta[0] = displacements[2 * elements[i].nodesIds[0] + 0];
      delta[1] = displacements[2 * elements[i].nodesIds[0] + 1];
      delta[2] = displacements[2 * elements[i].nodesIds[1] + 0];
      delta[3] = displacements[2 * elements[i].nodesIds[1] + 1];
      delta[4] = displacements[2 * elements[i].nodesIds[2] + 0];
      delta[5] = displacements[2 * elements[i].nodesIds[2] + 1];
    } else if (DIM == 3) {
      delta[0] = displacements[3 * elements[i].nodesIds[0] + 0];
      delta[1] = displacements[3 * elements[i].nodesIds[0] + 1];
      delta[2] = displacements[3 * elements[i].nodesIds[0] + 2];
      delta[3] = displacements[3 * elements[i].nodesIds[1] + 0];
      delta[4] = displacements[3 * elements[i].nodesIds[1] + 1];
      delta[5] = displacements[3 * elements[i].nodesIds[1] + 2];
      delta[6] = displacements[3 * elements[i].nodesIds[2] + 0];
      delta[7] = displacements[3 * elements[i].nodesIds[2] + 1];
      delta[8] = displacements[3 * elements[i].nodesIds[2] + 2];
      delta[9] = displacements[3 * elements[i].nodesIds[3] + 0];
      delta[10] = displacements[3 * elements[i].nodesIds[3] + 1];
      delta[11] = displacements[3 * elements[i].nodesIds[3] + 2];
    }

    for (int k = 0; k < 3 * (DIM - 1); ++k)
      for (int j = 0; j < 6 * (DIM - 1); ++j)
        B(k, j) = all_B[j + 6 * (DIM - 1) * k + 3 * (DIM - 1) * 6 * (DIM - 1) * i];

    //DeformationVector = it->B.Product(delta);
    DeformationVector = B.Product(delta);
    StressVector = D.Product(DeformationVector);


//    double sigma = sqrt(StressVector[0] * StressVector[0] - StressVector[0]
//        * StressVector[1] + StressVector[1] * StressVector[1] + 3.0 * StressVector[2] * StressVector[2]);
//    sigma_mises.push_back(sigma);

//    double epsilon = sqrt(DeformationVector[0] * DeformationVector[0] - DeformationVector[0]
//        * DeformationVector[1] + DeformationVector[1] * DeformationVector[1] + 3.0 * DeformationVector[2] * DeformationVector[2]);
//    epsilon_mises.push_back(epsilon);

    Deformation.push_back(DeformationVector);
    Stress.push_back(StressVector);

    //DeformationVector.Show();
    //StressVector.Show();
    //cout << endl;
  }
}

bool CheckPointInside(float v1, float v2, float x1, float y1, float x2, float y2, float x3, float y3) {
  float res1 = (y1 - y2) * v1 + (x2 - x1) * v2 + (x1 * y2 - y1 * x2);
  float res2 = (y2 - y3) * v1 + (x3 - x2) * v2 + (x2 * y3 - y2 * x3);
  float res3 = (y3 - y1) * v1 + (x1 - x3) * v2 + (x3 * y1 - y3 * x1);

  double eps = 1e-15;
  if ((res1 < 0.0 || std::abs(res1) < eps) &&
      (res2 < 0.0 || std::abs(res2) < eps) &&
      (res3 < 0.0 || std::abs(res3) < eps))
    return true;
  if ((res1 > 0.0 || std::abs(res1) < eps) &&
      (res2 > 0.0 || std::abs(res2) < eps) &&
      (res3 > 0.0 || std::abs(res3) < eps))
    return true;

  //    if ((res1 < 0.0) &&
  //            (res2 < 0.0) &&
  //            (res3 < 0.0))
  //        return true;
  //    if ((res1 > 0.0) &&
  //            (res2 > 0.0) &&
  //            (res3 > 0.0))
  //        return true;

  return false;
}

bool CheckPointInside2(float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {
  float b1 = (x2 * y3 - y2 * x3) / (x2 - x3);
  float b2 = (x1 * y3 - y1 * x3) / (x1 - x3);
  float b3 = (x1 * y2 - y1 * x2) / (x1 - x2);
  float k1 = (y2 - b1) / x2;
  float k2 = (y1 - b1) / x1;
  float k3 = (y2 - b1) / x2;

  float res1_1 = y0 - x0 * k1 - b1;
  float res1_2 = y1 - x1 * k1 - b1;

  float res2_1 = y0 - x0 * k2 - b2;
  float res2_2 = y2 - x2 * k2 - b2;

  float res3_1 = y0 - x0 * k3 - b3;
  float res3_2 = y3 - x3 * k3 - b3;

  bool p1 = res1_1 > 0 && res1_2 > 0;
  bool n1 = res1_1 < 0 && res1_2 < 0;

  bool p2 = res2_1 > 0 && res2_2 > 0;
  bool n2 = res2_1 < 0 && res2_2 < 0;

  bool p3 = res3_1 > 0 && res3_2 > 0;
  bool n3 = res3_1 < 0 && res3_2 < 0;

  if (p1 && p2 && p3)
    return true;
  if (p1 && p2 && n3)
    return true;
  if (p1 && n2 && p3)
    return true;
  if (p1 && n2 && n3)
    return true;

  if (n1 && p2 && p3)
    return true;
  if (n1 && p2 && n3)
    return true;
  if (n1 && n2 && p3)
    return true;
  if (n1 && n2 && n3)
    return true;


  return false;
}

//2D only
void CalculateStressAlongAxis(std::vector<float> &StressComponents,
                              std::string axe,
                              std::string stress_component,
                              float fixed_value,
                              float a,
                              float b,
                              std::vector<MyArray> &Stress,
                              std::vector<MyArray> &nodes,
                              std::vector<Element> elements) {
  CheckRunTime(__func__)

  int component_id = 0;
  if (stress_component == "xx") {
    component_id = 0;
  } else if (stress_component == "yy") {
    component_id = 1;
  } else if (stress_component == "xy") {
    component_id = 2;
  }

  MyArray range(100, a, b);
  for (int i = 0; i < range.get_size(); ++i) {
    float x = (axe == "x") ? range[i] : fixed_value;
    float y = (axe == "y") ? range[i] : fixed_value;
    for (int j = 0; j < elements.size(); ++j) {
      float x1 = nodes[0][elements[j].nodesIds[0]], y1 = nodes[1][elements[j].nodesIds[0]];
      float x2 = nodes[0][elements[j].nodesIds[1]], y2 = nodes[1][elements[j].nodesIds[1]];
      float x3 = nodes[0][elements[j].nodesIds[2]], y3 = nodes[1][elements[j].nodesIds[2]];
      if (CheckPointInside(x, y, x1, y1, x2, y2, x3, y3)) {
        StressComponents.push_back(range[i]);
        StressComponents.push_back(Stress[j][component_id]);
        //Stress[j][component_id] = 1e+6;
        break;
      }
    }
  }
}

float Interpolate(float x, float y, float x1, float y1, float x2, float y2, float x3, float y3, float f1, float f2, float f3) {

  float a1 = x2 * y3 - x3 * y2;
  float b1 = y2 - y3;
  float c1 = x3 - x2;

  float a2 = x3 * y1 - x1 * y3;
  float b2 = y3 - y1;
  float c2 = x1 - x3;

  float a3 = x1 * y2 - x2 * y1;
  float b3 = y1 - y2;
  float c3 = x2 - x1;

  float dlt = std::abs((x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3));

  float N1 = (a1 + x * b1 + y * c1) / dlt;
  float N2 = (a2 + x * b2 + y * c2) / dlt;
  float N3 = (a3 + x * b3 + y * c3) / dlt;

  return f1 * N1 + f2 * N2 + f3 * N3;
}

//2D only
void CalculateStressAlongAxisSmooth(std::vector<float> &StressComponentsSmooth,
                                    std::string axe,
                                    float fixed_value,
                                    float a,
                                    float b,
                                    MyArray StressSmooth,
                                    std::vector<MyArray> &nodes,
                                    std::vector<Element> elements) {
  CheckRunTime(__func__)

  MyArray range(100, a, b);
  for (int i = 0; i < range.get_size(); ++i) {
    float x = (axe == "x") ? range[i] : fixed_value;
    float y = (axe == "y") ? range[i] : fixed_value;
    for (int j = 0; j < elements.size(); ++j) {
      float x1 = nodes[0][elements[j].nodesIds[0]], y1 = nodes[1][elements[j].nodesIds[0]];
      float x2 = nodes[0][elements[j].nodesIds[1]], y2 = nodes[1][elements[j].nodesIds[1]];
      float x3 = nodes[0][elements[j].nodesIds[2]], y3 = nodes[1][elements[j].nodesIds[2]];
      if (CheckPointInside(x, y, x1, y1, x2, y2, x3, y3)) {
        StressComponentsSmooth.push_back(range[i]);
        float f1 = StressSmooth[elements[j].nodesIds[0]];
        float f2 = StressSmooth[elements[j].nodesIds[1]];
        float f3 = StressSmooth[elements[j].nodesIds[2]];
        StressComponentsSmooth.push_back(Interpolate(x, y, x1, y1, x2, y2, x3, y3, f1, f2, f3));
        break;
      }
    }
  }
}

void CalculateMisesAlongLineMises(std::vector<float> &MisesComponents,
                                  float k, float m,
                                  float a, float b,
                                  std::vector<float> sigma_mises,
                                  std::vector<MyArray> &nodes,
                                  std::vector<Element> elements) {
  CheckRunTime(__func__)

  MyArray range(100, a, b);
  for (int i = 0; i < range.get_size(); ++i) {
    float x = range[i];
    float y = k * x + m;
    for (int j = 0; j < elements.size(); ++j) {
      float x1 = nodes[0][elements[j].nodesIds[0]], y1 = nodes[1][elements[j].nodesIds[0]];
      float x2 = nodes[0][elements[j].nodesIds[1]], y2 = nodes[1][elements[j].nodesIds[1]];
      float x3 = nodes[0][elements[j].nodesIds[2]], y3 = nodes[1][elements[j].nodesIds[2]];
      if (CheckPointInside(x, y, x1, y1, x2, y2, x3, y3)) {
        MisesComponents.push_back(range[i]);
        MisesComponents.push_back(sigma_mises[j]);
        break;
      }
    }
  }
}

void CalculateNodeAdjElem(FEMdataKeeper FEMdata, std::unordered_map <int, std::vector<int>> &a) {
  for (int n = 0; n < FEMdata.elements.size(); ++n) {
    a[FEMdata.elements[n].nodesIds[0]].push_back(n);
    a[FEMdata.elements[n].nodesIds[1]].push_back(n);
    a[FEMdata.elements[n].nodesIds[2]].push_back(n);
  }
}

void ApplyConstraints_EbE(FEMdataKeeper &FEMdata) {
  CheckRunTime(__func__)
  std::vector<int> indicesToConstraint;
  for (std::vector<Constraint>::const_iterator it = FEMdata.constraints.begin(); it != FEMdata.constraints.end(); ++it) {
    if (it->type & Constraint::UX) {
      indicesToConstraint.push_back(FEMdata.DIM * it->node + 0);
    }
    if (it->type & Constraint::UY) {
      indicesToConstraint.push_back(FEMdata.DIM * it->node + 1);
    }
    if (it->type & Constraint::UZ) {
      indicesToConstraint.push_back(FEMdata.DIM * it->node + 2);
    }
  }

  //CUDA
  //std::cout << "CONSTRAINTS\n";
  FEMdata.CudaIndicesToConstraints.Resize(indicesToConstraint.size());
  for (int i = 0; i < indicesToConstraint.size(); ++i) {
    FEMdata.CudaIndicesToConstraints[i] = indicesToConstraint[i];
    //std::cout << FEMdata.CudaIndicesToConstraints[i] << " ";
  }
  FEMdata.CudaIndicesToConstraintsCount = indicesToConstraint.size();
  //

//  for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
//    for (int c_id = 0; c_id < indicesToConstraint.size(); c_id++) {
//      for (int i = 0; i < 3; ++i) {
//        for (int j = 0; j < 3; ++j) {
//          for (int ilocal = 0; ilocal < 2; ++ilocal) {
//            for (int jlocal = 0; jlocal < 2; ++jlocal) {
//              if (2 * it->nodesIds[i] + ilocal == indicesToConstraint[c_id] ||
//                  2 * it->nodesIds[j] + jlocal == indicesToConstraint[c_id]) {
//                if (2 * it->nodesIds[i] + ilocal != 2 * it->nodesIds[j] + jlocal) {
//                  it->Klocal(2 * i + ilocal, 2 * j + jlocal) = 0.0;
//                }
//              }
//            }
//          }
//        }
//      }
//    }
//  }
}

void PCG_EbE2(FEMdataKeeper &FEMdata, std::unordered_map <int, std::vector<int>> &node2adj_elem, float eps) {
  float gamma0 = 0.0f;
  float gamma = 0.0f;
  float gamma_new = 0.0f;

  int enum_it = 0;
  for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
    // (0a)
    it->r = it->Flocal;

    // (0b)
    for (int local_node = 0; local_node < 3; ++local_node) {
      std::vector<int> adjElems = node2adj_elem[it->nodesIds[local_node]];
      for (int i = 0; i < adjElems.size(); ++i) {
        int tmp = FEMdata.elements[adjElems[i]].Global2LocalNode(it->nodesIds[local_node]);
        it->m[local_node * 2 + 0] += float(FEMdata.elements[adjElems[i]].Klocal(tmp * 2 + 0, tmp * 2 + 0));
        it->m[local_node * 2 + 1] += float(FEMdata.elements[adjElems[i]].Klocal(tmp * 2 + 1, tmp * 2 + 1));
      }
    }

    // (0c)
    for (int local_node = 0; local_node < 3; ++local_node) {
      it->z[local_node * 2 + 0] = it->r[local_node * 2 + 0] / it->m[local_node * 2 + 0];
      it->z[local_node * 2 + 1] = it->r[local_node * 2 + 1] / it->m[local_node * 2 + 1];
    }

    ++enum_it;

  }
}

void GenerateMask(FEMdataKeeper FEMdata, MyArray &mask) {
  CheckRunTime(__func__)
  int DIM = FEMdata.DIM;
  for (int eIdx = 0; eIdx < FEMdata.elementsCount; ++eIdx) {
    mask[3 * DIM * eIdx + 0] = FEMdata.elements[eIdx].nodesIds[0] * DIM + 0;
    mask[3 * DIM * eIdx + 1] = FEMdata.elements[eIdx].nodesIds[0] * DIM + 1;
    mask[3 * DIM * eIdx + 2] = FEMdata.elements[eIdx].nodesIds[1] * DIM + 0;
    mask[3 * DIM * eIdx + 3] = FEMdata.elements[eIdx].nodesIds[1] * DIM + 1;
    mask[3 * DIM * eIdx + 4] = FEMdata.elements[eIdx].nodesIds[2] * DIM + 0;
    mask[3 * DIM * eIdx + 5] = FEMdata.elements[eIdx].nodesIds[2] * DIM + 1;
  }
}

void solve_diag_vec(FEMdataKeeper &FEMdata, MyArray &res, bool doAssemblyRes) {
  CheckRunTime(__func__)
  int DIM = FEMdata.DIM;
  int n_elems  = FEMdata.elementsCount;
  int n_gl_dofs = FEMdata.nodesCount * DIM;

  int grid_size = 3 * DIM * n_elems;
  MyArray r(grid_size), diag(grid_size), mask(grid_size),
      diag_assemblied(n_gl_dofs), r_assemblied(n_gl_dofs);

  //GenerateMask(FEMdata, mask);

  MyArray n_adjelem(n_gl_dofs);

  SparseMatrixCOO Q (grid_size);          // Reduction mask Q

  for (int eIdx = 0; eIdx < n_elems; ++eIdx) {
    Element elem = FEMdata.elements[eIdx];
    Q.write_value(3 * DIM * eIdx + 0 * DIM + 0, elem.nodesIds[0] * DIM + 0, 1);
    Q.write_value(3 * DIM * eIdx + 0 * DIM + 1, elem.nodesIds[0] * DIM + 1, 1);
    Q.write_value(3 * DIM * eIdx + 1 * DIM + 0, elem.nodesIds[1] * DIM + 0, 1);
    Q.write_value(3 * DIM * eIdx + 1 * DIM + 1, elem.nodesIds[1] * DIM + 1, 1);
    Q.write_value(3 * DIM * eIdx + 2 * DIM + 0, elem.nodesIds[2] * DIM + 0, 1);
    Q.write_value(3 * DIM * eIdx + 2 * DIM + 1, elem.nodesIds[2] * DIM + 1, 1);

    n_adjelem[elem.nodesIds[0] * DIM + 0] += 1.0f;
    n_adjelem[elem.nodesIds[0] * DIM + 1] += 1.0f;
    n_adjelem[elem.nodesIds[1] * DIM + 0] += 1.0f;
    n_adjelem[elem.nodesIds[1] * DIM + 1] += 1.0f;
    n_adjelem[elem.nodesIds[2] * DIM + 0] += 1.0f;
    n_adjelem[elem.nodesIds[2] * DIM + 1] += 1.0f;

    diag[3 * DIM * eIdx + 0] = elem.Alocal(0, 0);
    diag[3 * DIM * eIdx + 1] = elem.Alocal(1, 1);
    diag[3 * DIM * eIdx + 2] = elem.Alocal(2, 2);
    diag[3 * DIM * eIdx + 3] = elem.Alocal(3, 3);
    diag[3 * DIM * eIdx + 4] = elem.Alocal(4, 4);
    diag[3 * DIM * eIdx + 5] = elem.Alocal(5, 5);

    r[3 * DIM * eIdx + 0] = elem.blocal[0];
    r[3 * DIM * eIdx + 1] = elem.blocal[1];
    r[3 * DIM * eIdx + 2] = elem.blocal[2];
    r[3 * DIM * eIdx + 3] = elem.blocal[3];
    r[3 * DIM * eIdx + 4] = elem.blocal[4];
    r[3 * DIM * eIdx + 5] = elem.blocal[5];
  }

  // на листочке распиши!
  diag_assemblied = Q.MultiplyTransposedByVector(diag, n_gl_dofs);
  r_assemblied = Q.MultiplyTransposedByVector(r, n_gl_dofs);

  MyArray x(n_gl_dofs);
  x = r_assemblied.divideByElementwise(diag_assemblied);

  if (doAssemblyRes) {
    res = x;
  } else {
    res = Q.MultiplyByVector(x, grid_size);
  }
}

void PCG_EbE_vec(FEMdataKeeper &FEMdata, MyArray &res, bool doAssemblyRes, float eps, bool PRINT_DEBUG_INFO) {
  CheckRunTime(__func__)
  int DIM = FEMdata.DIM;
  int n_elems  = FEMdata.elementsCount;
  int n_gl_dofs = FEMdata.nodesCount * DIM;

  int grid_size = 3 * DIM * n_elems;
  MyArray r(grid_size), diag(grid_size), m(grid_size),
      z(grid_size), s(grid_size), p(grid_size),
      u(grid_size), x(grid_size), mask(grid_size);

  GenerateMask(FEMdata, mask);

  MyArray n_adjelem(n_gl_dofs);

  SparseMatrixCOO Q (grid_size);          // Reduction mask Q
  SparseMatrixCOO Ae(FEMdata.nonzeroMatrixNumCount);    // Block diagonal matrix consisting of A^(e)
  for (int eIdx = 0; eIdx < n_elems; ++eIdx) {
    Element elem = FEMdata.elements[eIdx];
    Q.write_value(3 * DIM * eIdx + 0 * DIM + 0, elem.nodesIds[0] * DIM + 0, 1);
    Q.write_value(3 * DIM * eIdx + 0 * DIM + 1, elem.nodesIds[0] * DIM + 1, 1);
    Q.write_value(3 * DIM * eIdx + 1 * DIM + 0, elem.nodesIds[1] * DIM + 0, 1);
    Q.write_value(3 * DIM * eIdx + 1 * DIM + 1, elem.nodesIds[1] * DIM + 1, 1);
    Q.write_value(3 * DIM * eIdx + 2 * DIM + 0, elem.nodesIds[2] * DIM + 0, 1);
    Q.write_value(3 * DIM * eIdx + 2 * DIM + 1, elem.nodesIds[2] * DIM + 1, 1);

    if (doAssemblyRes) {
      n_adjelem[elem.nodesIds[0] * DIM + 0] += 1.0f;
      n_adjelem[elem.nodesIds[0] * DIM + 1] += 1.0f;
      n_adjelem[elem.nodesIds[1] * DIM + 0] += 1.0f;
      n_adjelem[elem.nodesIds[1] * DIM + 1] += 1.0f;
      n_adjelem[elem.nodesIds[2] * DIM + 0] += 1.0f;
      n_adjelem[elem.nodesIds[2] * DIM + 1] += 1.0f;
    }

    diag[3 * DIM * eIdx + 0] = elem.Alocal(0, 0);
    diag[3 * DIM * eIdx + 1] = elem.Alocal(1, 1);
    diag[3 * DIM * eIdx + 2] = elem.Alocal(2, 2);
    diag[3 * DIM * eIdx + 3] = elem.Alocal(3, 3);
    diag[3 * DIM * eIdx + 4] = elem.Alocal(4, 4);
    diag[3 * DIM * eIdx + 5] = elem.Alocal(5, 5);

    for (int xlocal = 0; xlocal < 3 * DIM; ++xlocal) {
      for (int ylocal = 0; ylocal < 3 * DIM; ++ylocal) {
        Ae.write_value(3 * DIM * eIdx + xlocal, 3 * DIM * eIdx + ylocal, elem.Alocal(xlocal, ylocal));
      }
    }

    // (0a)
    // Initialize vector r^(e)
    r[3 * DIM * eIdx + 0] = elem.blocal[0];
    r[3 * DIM * eIdx + 1] = elem.blocal[1];
    r[3 * DIM * eIdx + 2] = elem.blocal[2];
    r[3 * DIM * eIdx + 3] = elem.blocal[3];
    r[3 * DIM * eIdx + 4] = elem.blocal[4];
    r[3 * DIM * eIdx + 5] = elem.blocal[5];
  }

  // (0b)

  //m = Q.MultiplyByVector(Q.MultiplyTransposedByVector(diag, n_gl_dofs), grid_size);
  gpuReductionWithMaskAndTransform(diag.get_data(), mask.get_data(), grid_size, m.get_data(), n_gl_dofs);

  //    std::cout << "CPU\n";
  //    diag.Show();
  //    mask.Show();
  //m.Show();
  //exit(-1);

  // (0c)
  z = r.divideByElementwise(m);

  //    diag.Show();
  //    r.Show();

  //m.Show();
  //z.Show();

  // (0d)
  //s = Q.MultiplyByVector(Q.MultiplyTransposedByVector(z, n_gl_dofs), grid_size);
  gpuReductionWithMaskAndTransform(z.get_data(), mask.get_data(), grid_size, s.get_data(), n_gl_dofs);

  //float gamma0 = r.dot_product(s);
  float gamma0 = gpuDotProduct(r.get_data(), s.get_data(), grid_size);
  float gamma = gamma0;
  float gamma_new = 0.0f;

  // (0e)
  p = s;
  if (PRINT_DEBUG_INFO) {
    //std::cout.precision(16);
    std::cout << "gamma0\t\t\t= " << gamma0 << std::endl << std::endl;
  }
  int n_iter = 0;
  do {
    ++n_iter;

    // (1a)

    u = Ae.MultiplyByVector(p, grid_size);
    //        gpuCalculateKlocal(FEMdata.CudaElements.get_data(), FEMdata.elementsCount,
    //                           FEMdata.nodesX.get_data(), FEMdata.nodesY.get_data(), FEMdata.nodesCount,
    //                           p.get_data(), u.get_data(), FEMdata.D.get_data(),
    //                           FEMdata.CudaIndicesToConstraints.get_data(), FEMdata.CudaIndicesToConstraintsCount);

    // (1b)
    float sumElem = p.dot_product(u);
    //float sumElem = gpuDotProduct(p.get_data(), u.get_data(), grid_size);
    // (1c,d)
    float alpha = gamma / sumElem;

    // (2a)
    x.add_weighted(p, 1.0f, alpha);
    //gpuAddWeighted(x.get_data(), p.get_data(), 1.0f, alpha, grid_size);



    // (2b)
    r.add_weighted(u, 1.0f, -1.0f * alpha);
    //gpuAddWeighted(r.get_data(), u.get_data(), 1.0f, -1.0f * alpha, grid_size);

    // (3)
    z = r.divideByElementwise(m);

    // (4)
    //s = Q.MultiplyByVector(Q.MultiplyTransposedByVector(z, n_gl_dofs), grid_size);
    gpuReductionWithMaskAndTransform(z.get_data(), mask.get_data(), grid_size, s.get_data(), n_gl_dofs);


    gamma_new = r.dot_product(s);
    //gamma_new = gpuDotProduct(r.get_data(), s.get_data(), grid_size);

    if (PRINT_DEBUG_INFO) {
      // Verbose
      // -----------------------------------------------------------------------
      std::cout << "Iteration #" << n_iter << std::endl;
      std::cout << "alpha (gamma / sumElem)\t= " << alpha << std::endl;
      std::cout << "alpha numerator (gamma)\t= " << gamma << std::endl;
      std::cout << "alpha denominator\t= " << sumElem << std::endl;
      // -----------------------------------------------------------------------
    }
    // (5)
    if (gamma_new < eps * gamma0)
      break;

    // (6)
    p.add_weighted(s, gamma_new / gamma, 1.0f);
    //gpuAddWeighted(p.get_data(), s.get_data(), gamma_new / gamma, 1.0f, grid_size);
    if (PRINT_DEBUG_INFO) {
      // Verbose
      // -----------------------------------------------------------------------
      std::cout << "beta\t\t\t= " << gamma_new / gamma << std::endl;
      std::cout << "gamma_new\t\t= " << gamma_new << std::endl;
      std::cout << std::endl;
      // -----------------------------------------------------------------------
    }
    gamma = gamma_new;

  } while (1);
  if (doAssemblyRes) {
    res = Q.MultiplyTransposedByVector(x, n_gl_dofs); //TODO: JUST REDUCTION
    //gpuReductionWithMask(x.get_data(), mask.get_data(), grid_size, res.get_data(), n_gl_dofs);
    res = res.divideByElementwise(n_adjelem);
  } else {
    res = x;
  }

}

void PCG_EbE_vec_cpu(FEMdataKeeper &FEMdata, MyArray &res, bool doAssemblyRes, float eps) {
  CheckRunTime(__func__)
  int DIM = FEMdata.DIM;
  int n_elems  = FEMdata.elementsCount;
  int n_gl_dofs = FEMdata.nodesCount * DIM;

  int grid_size = 3 * DIM * n_elems;
  MyArray r(grid_size), diag(grid_size), m(grid_size),
      z(grid_size), s(grid_size), p(grid_size),
      u(grid_size), x(grid_size);

  MyArray n_adjelem(n_gl_dofs);

  SparseMatrixCOO Q (grid_size);          // Reduction mask Q
  SparseMatrixCOO Ae(FEMdata.nonzeroMatrixNumCount);    // Block diagonal matrix consisting of A^(e)
  for (int eIdx = 0; eIdx < n_elems; ++eIdx) {
    Element elem = FEMdata.elements[eIdx];
    Q.write_value(3 * DIM * eIdx + 0 * DIM + 0, elem.nodesIds[0] * DIM + 0, 1);
    Q.write_value(3 * DIM * eIdx + 0 * DIM + 1, elem.nodesIds[0] * DIM + 1, 1);
    Q.write_value(3 * DIM * eIdx + 1 * DIM + 0, elem.nodesIds[1] * DIM + 0, 1);
    Q.write_value(3 * DIM * eIdx + 1 * DIM + 1, elem.nodesIds[1] * DIM + 1, 1);
    Q.write_value(3 * DIM * eIdx + 2 * DIM + 0, elem.nodesIds[2] * DIM + 0, 1);
    Q.write_value(3 * DIM * eIdx + 2 * DIM + 1, elem.nodesIds[2] * DIM + 1, 1);

    n_adjelem[elem.nodesIds[0] * DIM + 0] += 1.0f;
    n_adjelem[elem.nodesIds[0] * DIM + 1] += 1.0f;
    n_adjelem[elem.nodesIds[1] * DIM + 0] += 1.0f;
    n_adjelem[elem.nodesIds[1] * DIM + 1] += 1.0f;
    n_adjelem[elem.nodesIds[2] * DIM + 0] += 1.0f;
    n_adjelem[elem.nodesIds[2] * DIM + 1] += 1.0f;

    diag[3 * DIM * eIdx + 0] = elem.Alocal(0, 0);
    diag[3 * DIM * eIdx + 1] = elem.Alocal(1, 1);
    diag[3 * DIM * eIdx + 2] = elem.Alocal(2, 2);
    diag[3 * DIM * eIdx + 3] = elem.Alocal(3, 3);
    diag[3 * DIM * eIdx + 4] = elem.Alocal(4, 4);
    diag[3 * DIM * eIdx + 5] = elem.Alocal(5, 5);

    for (int xlocal = 0; xlocal < 3 * DIM; ++xlocal) {
      for (int ylocal = 0; ylocal < 3 * DIM; ++ylocal) {
        Ae.write_value(3 * DIM * eIdx + xlocal, 3 * DIM * eIdx + ylocal, elem.Alocal(xlocal, ylocal));
      }
    }

    // (0a)
    // Initialize vector r^(e)
    r[3 * DIM * eIdx + 0] = elem.blocal[0];
    r[3 * DIM * eIdx + 1] = elem.blocal[1];
    r[3 * DIM * eIdx + 2] = elem.blocal[2];
    r[3 * DIM * eIdx + 3] = elem.blocal[3];
    r[3 * DIM * eIdx + 4] = elem.blocal[4];
    r[3 * DIM * eIdx + 5] = elem.blocal[5];
  }

  // (0b)
  m = Q.MultiplyByVector(Q.MultiplyTransposedByVector(diag, n_gl_dofs), grid_size);

  // (0c)
  z = r.divideByElementwise(m);

  // (0d)
  s = Q.MultiplyByVector(Q.MultiplyTransposedByVector(z, n_gl_dofs), grid_size);
  float gamma0 = r.dot_product(s);
  float gamma = gamma0;
  float gamma_new = 0.0f;

  // (0e)
  p = s;

  std::cout.precision(16);
  std::cout << "gamma0\t\t\t= " << gamma0 << std::endl << std::endl;

  int n_iter = 0;
  do {
    ++n_iter;

    // (1a)
    u = Ae.MultiplyByVector(p, grid_size);
    // (1b)
    float sumElem = p.dot_product(u);
    // (1c,d)
    float alpha = gamma / sumElem;

    // (2a)
    x.add_weighted(p, 1.0f, alpha);
    // (2b)
    r.add_weighted(u, 1.0f, -1.0f * alpha);

    // (3)
    z = r.divideByElementwise(m);

    // (4)
    s = Q.MultiplyByVector(Q.MultiplyTransposedByVector(z, n_gl_dofs), grid_size);
    gamma_new = r.dot_product(s);

    // Verbose
    // -----------------------------------------------------------------------
    std::cout << "Iteration #" << n_iter << std::endl;
    std::cout << "alpha (gamma / sumElem)\t= " << alpha << std::endl;
    std::cout << "alpha numerator (gamma)\t= " << gamma << std::endl;
    std::cout << "alpha denominator\t= " << sumElem << std::endl;
    // -----------------------------------------------------------------------

    // (5)
    if (gamma_new < eps * gamma0)
      break;

    // (6)
    p.add_weighted(s, gamma_new / gamma, 1.0f);

    // Verbose
    // -----------------------------------------------------------------------
    std::cout << "beta\t\t\t= " << gamma_new / gamma << std::endl;
    std::cout << "gamma_new\t\t= " << gamma_new << std::endl;
    std::cout << std::endl;
    // -----------------------------------------------------------------------

    gamma = gamma_new;

  } while (1);

  if (doAssemblyRes) {
    res = Q.MultiplyTransposedByVector(x, n_gl_dofs);
    res = res.divideByElementwise(n_adjelem);
  } else
    res = x;

}

// Note that instead of Klocal is Alocal
void PCG_EbE(FEMdataKeeper &FEMdata, std::unordered_map <int, std::vector<int>> &node2adj_elem, float eps, bool PRINT_DEBUG_INFO) {
  // see article "A distributed memory parallel element-by-element scheme based on Jacobi-conditioned conjugate gradient for 3D finite element analysis"
  // by Yaoru Liu, Weiyuan Zhou, Qiang Yang
  CheckRunTime(__func__)

  float gamma0 = 0.0f;
  float gamma = 0.0f;
  float gamma_new = 0.0f;

  int enum_it = 0;

  for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
    // (0a)
    it->r = it->blocal;

    // (0b)
    for (int local_node = 0; local_node < 3; ++local_node) {
      std::vector<int> adjElems = node2adj_elem[it->nodesIds[local_node]];
      for (int i = 0; i < adjElems.size(); ++i) {
        int tmp = FEMdata.elements[adjElems[i]].Global2LocalNode(it->nodesIds[local_node]);
        it->m[local_node * 2 + 0] += float(FEMdata.elements[adjElems[i]].Alocal(tmp * 2 + 0, tmp * 2 + 0));
        it->m[local_node * 2 + 1] += float(FEMdata.elements[adjElems[i]].Alocal(tmp * 2 + 1, tmp * 2 + 1));
      }
    }

    // (0c)
    for (int local_node = 0; local_node < 3; ++local_node) {
      it->z[local_node * 2 + 0] = it->r[local_node * 2 + 0] / it->m[local_node * 2 + 0];
      it->z[local_node * 2 + 1] = it->r[local_node * 2 + 1] / it->m[local_node * 2 + 1];
    }

    ++enum_it;

  }

  enum_it = 0;
  for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
    // (0d)
    for (int local_node = 0; local_node < 3; ++local_node) {
      it->s[local_node * 2 + 0] = 0.0f;
      it->s[local_node * 2 + 1] = 0.0f;
      std::vector<int> adjElems = node2adj_elem[it->nodesIds[local_node]];
      for (int i = 0; i < adjElems.size(); ++i) {
        int tmp = FEMdata.elements[adjElems[i]].Global2LocalNode(it->nodesIds[local_node]);
        it->s[local_node * 2 + 0] += FEMdata.elements[adjElems[i]].z[tmp * 2 + 0];
        it->s[local_node * 2 + 1] += FEMdata.elements[adjElems[i]].z[tmp * 2 + 1];
      }
    }
    gamma0 += it->r.dot_product(it->s);

    // (0e)
    it->p = it->s;

    ++enum_it;
  }

  gamma = gamma0;

  if (PRINT_DEBUG_INFO) {
    //std::cout.precision(16);
    std::cout << "gamma0\t\t\t= " << gamma0 << std::endl << std::endl;
  }

  int n_iter = 0;
  do {
    ++n_iter;
    int enum_it = 0;
    if (PRINT_DEBUG_INFO) {
      std::cout << "Iteration #" << n_iter << std::endl;
    }
    float sumElem = 0.0f;
    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
      // (1a-d)
      it->u = it->Alocal.Product(it->p);
      float tmp = it->p.dot_product(it->u);
      sumElem += tmp;
    }

    float alpha = gamma / sumElem;

    if (PRINT_DEBUG_INFO) {
      std::cout << "alpha (gamma / sumElem)\t= " << alpha << std::endl;
      std::cout << "alpha numerator (gamma)\t= " << gamma << std::endl;
      std::cout << "alpha denominator\t= " << sumElem << std::endl;
    }

    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
      // (2a)
      it->res.add_weighted(it->p, 1.0f, alpha);
      // (2b)
      it->r.add_weighted(it->u, 1.0f, -1.0f * alpha);

      // (3)
      for (int local_node = 0; local_node < 3; ++local_node) {
        it->z[local_node * 2 + 0] = it->r[local_node * 2 + 0] / it->m[local_node * 2 + 0];
        it->z[local_node * 2 + 1] = it->r[local_node * 2 + 1] / it->m[local_node * 2 + 1];
      }
    }

    gamma_new = 0.0f;
    // (4 a-e)
    enum_it = 0;
    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
      for (int local_node = 0; local_node < 3; ++local_node) {
        it->s[local_node * 2 + 0] = 0.0f;
        it->s[local_node * 2 + 1] = 0.0f;
        std::vector<int> adjElems = node2adj_elem[it->nodesIds[local_node]];
        for (int i = 0; i < adjElems.size(); ++i) {
          int tmp = FEMdata.elements[adjElems[i]].Global2LocalNode(it->nodesIds[local_node]);
          it->s[local_node * 2 + 0] += FEMdata.elements[adjElems[i]].z[tmp * 2 + 0];
          it->s[local_node * 2 + 1] += FEMdata.elements[adjElems[i]].z[tmp * 2 + 1];
        }
      }

      gamma_new += it->r.dot_product(it->s);

      ++enum_it;
    }

    // (5)
    if (gamma_new < eps * gamma0)
      break;

    // (6)
    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
      it->p.add_weighted(it->s, gamma_new / gamma, 1.0f);
    }

    if (PRINT_DEBUG_INFO) {
      std::cout << "beta\t\t\t= " << gamma_new / gamma << std::endl;
      std::cout << "gamma_new\t\t= " << gamma_new << std::endl;
    }

    gamma = gamma_new;

  } while (1);
}

void writeSnapshot(float t, int num_receivers, int grid_size, int n_gl_dofs, FEMdataKeeper &FEMdata, gpuDataKeeper_DYN &gpu_data) {
  MyArray h_temp(n_gl_dofs);
  gpuReductionWithMask2(gpu_data.get_x(), gpu_data.get_mask(), grid_size, gpu_data.get_displ_global());
  gpuCopyDeviceToHost(gpu_data.get_displ_global(), h_temp.get_data(), n_gl_dofs);

  fstream out;
  out.open(FEMdata.res_dir + "/snapshots.csv", fstream::out | fstream::app);
  out << t;
  for (int i = 1; i < 2 * num_receivers + 1; i += 2)
    out << " " << h_temp[i] << " " << h_temp[i + 1];
  out << "\n";
  out.close();
}

void CalculateFEM(FEMdataKeeper &FEMdata, bool PRINT_DEBUG_INFO) {
  CheckRunTime(__func__)
  for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
    it->CalculateKlocal(FEMdata.D, FEMdata.nodes);
  }

  for (std::vector<BoundaryEdge>::iterator it = FEMdata.boundary.begin(); it != FEMdata.boundary.end(); ++it) {
    if (FEMdata.DIM == 2)
      FEMdata.elements[it->adj_elem1].CalculateFlocal2D(*it, FEMdata.nodes, 0.0f);
    else if (FEMdata.DIM == 3)
      FEMdata.elements[it->adj_elem1].CalculateFlocal3D(*it, FEMdata.nodes, 0.0f);
  }

  SparseMatrixCOO globalK = AssemblyStiffnessMatrix   (FEMdata);

  MyArray F = AssemblyF(FEMdata); // globalK * displacements = F
  //F.add(FEMdata.loads);
  ApplyLoads(F, FEMdata.loads);
  ApplyConstraints(globalK, F, FEMdata.constraints); //, FEMdata.nodesCount);

  globalK.SortIt();

  globalK.PCG_solve(F, FEMdata.displacements, 1e-10f);
}

void CalculateFEM_EbE(FEMdataKeeper &FEMdata, bool PRINT_DEBUG_INFO) {
  CheckRunTime(__func__)
  for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
    it->CalculateKlocal(FEMdata.D, FEMdata.nodes);
  }

  for (std::vector<BoundaryEdge>::iterator it = FEMdata.boundary.begin(); it != FEMdata.boundary.end(); ++it) {
    if (FEMdata.DIM == 2)
      FEMdata.elements[it->adj_elem1].CalculateFlocal2D(*it, FEMdata.nodes, 0.0f);
    else if (FEMdata.DIM == 3)
      FEMdata.elements[it->adj_elem1].CalculateFlocal3D(*it, FEMdata.nodes, 0.0f);
  }

  std::unordered_map <int, std::vector<int>> nodeAdjElem;
  CalculateNodeAdjElem(FEMdata, nodeAdjElem);

  ApplyConstraints_EbE(FEMdata);
  AssignLoadElement(FEMdata, nodeAdjElem);
  ApplyLoads_EbE(FEMdata);

  for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
    it->Alocal = it->Klocal;
  }
  for (std::vector<BoundaryEdge>::iterator it = FEMdata.boundary.begin(); it != FEMdata.boundary.end(); ++it) {
    FEMdata.elements[it->adj_elem1].blocal.add(FEMdata.elements[it->adj_elem1].Flocal);
  }

  PCG_EbE(FEMdata, nodeAdjElem, 1e-10f, PRINT_DEBUG_INFO);
  for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
    it->x = it->res;
  }
  AssemblyX(FEMdata, nodeAdjElem);
}

void CalculateFEM_EbE_vec(FEMdataKeeper &FEMdata, bool PRINT_DEBUG_INFO) {
  CheckRunTime(__func__)
  for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
    it->CalculateKlocal(FEMdata.D, FEMdata.nodes);
  }

  for (std::vector<BoundaryEdge>::iterator it = FEMdata.boundary.begin(); it != FEMdata.boundary.end(); ++it) {
    if (FEMdata.DIM == 2)
      FEMdata.elements[it->adj_elem1].CalculateFlocal2D(*it, FEMdata.nodes, 0.0f);
    else if (FEMdata.DIM == 3)
      FEMdata.elements[it->adj_elem1].CalculateFlocal3D(*it, FEMdata.nodes, 0.0f);
  }

  std::unordered_map <int, std::vector<int>> nodeAdjElem;
  CalculateNodeAdjElem(FEMdata, nodeAdjElem);

  AssignLoadElement(FEMdata, nodeAdjElem);
  ApplyLoads_EbE(FEMdata);

  int nonzero_matrix_numbers = 0;
  ApplyConstraints_EbE(FEMdata);
  for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
    nonzero_matrix_numbers += it->Klocal.CountNonzero();
  }
  FEMdata.nonzeroMatrixNumCount = nonzero_matrix_numbers;

  for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
    it->Alocal = it->Klocal;
  }
  for (std::vector<BoundaryEdge>::iterator it = FEMdata.boundary.begin(); it != FEMdata.boundary.end(); ++it) {
    FEMdata.elements[it->adj_elem1].blocal.add(FEMdata.elements[it->adj_elem1].Flocal);
  }
  PCG_EbE_vec(FEMdata, FEMdata.displacements, true, 1e-10f, PRINT_DEBUG_INFO);
}

void CalculateFEM_EbE_vec_GPU(FEMdataKeeper &FEMdata, bool PRINT_DEBUG_INFO) {
  CheckRunTime(__func__)
  //    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
  //        it->CalculateKlocal(FEMdata.D, FEMdata.nodes);
  //    }

  for (std::vector<BoundaryEdge>::iterator it = FEMdata.boundary.begin(); it != FEMdata.boundary.end(); ++it) {
    if (FEMdata.DIM == 2)
      FEMdata.elements[it->adj_elem1].CalculateFlocal2D(*it, FEMdata.nodes, 0.0f);
    else if (FEMdata.DIM == 3)
      FEMdata.elements[it->adj_elem1].CalculateFlocal3D(*it, FEMdata.nodes, 0.0f);
  }

  //    std::unordered_map <int, std::vector<int>> nodeAdjElem;
  //    CalculateNodeAdjElem(FEMdata, nodeAdjElem);

  //    AssignLoadElement(FEMdata, nodeAdjElem);
  //    ApplyLoads_EbE(FEMdata);

  //    int nonzero_matrix_numbers = 0;
  ApplyConstraints_EbE(FEMdata);
  //    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
  //        nonzero_matrix_numbers += it->Klocal.CountNonzero();
  //    }
  //    FEMdata.nonzeroMatrixNumCount = nonzero_matrix_numbers;

  //    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
  //        it->Alocal = it->Klocal;
  //    }
  //    for (std::vector<BoundaryEdge>::iterator it = FEMdata.boundary.begin(); it != FEMdata.boundary.end(); ++it) {
  //        FEMdata.elements[it->adj_elem1].blocal.add(FEMdata.elements[it->adj_elem1].Flocal);
  //    }
  //    PCG_EbE_vec(FEMdata, FEMdata.displacements, true, 1e-10f, PRINT_DEBUG_INFO);

  // Arseniy:
  gpuPCG_EbE_vec(FEMdata, FEMdata.displacements, true, 1e-4f, PRINT_DEBUG_INFO);

//  // Grisha 1:
//  bool doAssemblyRes = true;
//  gpuDataKeeper gpu_data(FEMdata.elementsCount, FEMdata.nodesCount, doAssemblyRes);
//  copyElements(FEMdata, gpu_data);
//  gpuCalculateKlocal2(gpu_data, FEMdata.elementsCount, FEMdata.nodesX.get_data(), FEMdata.nodesY.get_data(),
//                      FEMdata.nodesCount, FEMdata.D.get_data(), FEMdata.CudaIndicesToConstraints.get_data(),
//                      FEMdata.CudaIndicesToConstraintsCount);
//  gpuPCG_EbE_vec_gpudataArg(FEMdata, FEMdata.displacements, true, 1e-4f, PRINT_DEBUG_INFO, gpu_data);

//  // Grisha 2:
//  bool doAssemblyRes = true;
//  bool isLumped = false;
//  gpuDataKeeper_DYN gpu_data(FEMdata.elementsCount, FEMdata.nodesCount, doAssemblyRes, isLumped);
//  float beta2 = 0.5; float dt = 0.02;
//  gpu_data.set_SLAU_matrix_coefs(1.0f, 0.5f * beta2 * dt*dt);

//  copyElementsAndFlocals(FEMdata, gpu_data);
//  gpuCalculateKlocal2(gpu_data, FEMdata.elementsCount, FEMdata.nodesX.get_data(), FEMdata.nodesY.get_data(),
//                      FEMdata.nodesCount, FEMdata.D.get_data(), FEMdata.CudaIndicesToConstraints.get_data(),
//                      FEMdata.CudaIndicesToConstraintsCount);
//  gpuCalculateMlocal(gpu_data, FEMdata.elementsCount,
//                           FEMdata.nodesX.get_data(), FEMdata.nodesY.get_data(), FEMdata.nodesCount);

//  gpuPCG_EbE_vec_DYN(FEMdata, gpu_data, doAssemblyRes, 1e-4f, PRINT_DEBUG_INFO);

//  int n_gl_dofs = FEMdata.nodesCount * DIM;
//  gpuCopyDeviceToHost(gpu_data.get_temp_res(), FEMdata.displacements.get_data(), n_gl_dofs);
}

void gpuCalculateFEM_EbE_vec(FEMdataKeeper &FEMdata, bool PRINT_DEBUG_INFO) {
  CheckRunTime(__func__)
  //    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
  //        it->CalculateKlocal(FEMdata.D, FEMdata.nodes);
  //    }

  int num = 0;
  for (std::vector<BoundaryEdge>::iterator it = FEMdata.boundary.begin(); it != FEMdata.boundary.end(); ++it) {
    if (FEMdata.DIM == 2)
      FEMdata.elements[it->adj_elem1].CalculateFlocal2D(*it, FEMdata.nodes, FEMdata.pressure[num++]);
    else if (FEMdata.DIM == 3)
      FEMdata.elements[it->adj_elem1].CalculateFlocal3D(*it, FEMdata.nodes, FEMdata.pressure[num++]);
  }

  std::unordered_map <int, std::vector<int>> nodeAdjElem;
  CalculateNodeAdjElem(FEMdata, nodeAdjElem);

  AssignLoadElement(FEMdata, nodeAdjElem);
  ApplyLoads_EbE(FEMdata);


  gpuPCG_EbE_vec(FEMdata, FEMdata.displacements, true, 1e-10f, PRINT_DEBUG_INFO);
}

void gpuPCG_EbE_vec(FEMdataKeeper &FEMdata, MyArray &res, bool doAssemblyRes, float eps, bool PRINT_DEBUG_INFO) {
  CheckRunTime(__func__)
  int n_elems  = FEMdata.elementsCount;
  int n_gl_dofs = FEMdata.nodesCount * FEMdata.DIM;
  int grid_size = 6 * (FEMdata.DIM - 1) * n_elems;

  gpuDataKeeper gpu_data(FEMdata.DIM, FEMdata.elementsCount, FEMdata.nodesCount, doAssemblyRes);

  copyElementsAndFlocals(FEMdata, gpu_data);
  gpuGenerateMask(gpu_data, FEMdata.DIM, FEMdata.elementsCount);
  gpuCountNAdjElem(gpu_data, grid_size);

  gpuCalculateKlocal2(gpu_data, FEMdata);
  gpu_data.copyBmatrixToHost(FEMdata.all_B.get_data());
  // (0a)
  gpuCopyDeviceToDevice(gpu_data.get_Flocals(), gpu_data.get_r(), grid_size);
  // (0b)
  gpuReductionWithMaskAndTransform2(gpu_data.get_diag(), gpu_data.get_mask(), grid_size, gpu_data.get_m(), n_gl_dofs);
  // (0c)
  gpuDivideByElementwise(gpu_data.get_r(), gpu_data.get_m(), gpu_data.get_z(), grid_size);
  // (0d)
  gpuReductionWithMaskAndTransform2(gpu_data.get_z(), gpu_data.get_mask(), grid_size, gpu_data.get_s(), n_gl_dofs);
  float gamma0 = gpuDotProduct2(gpu_data.get_r(), gpu_data.get_s(), grid_size);
  float gamma = gamma0;
  float gamma_new = 0.0f;

  // (0e)
  gpuCopyDeviceToDevice(gpu_data.get_s(), gpu_data.get_p(), grid_size);

  if (PRINT_DEBUG_INFO) {
    //std::cout.precision(16);
    std::cout << "gamma0\t\t\t= " << gamma0 << std::endl << std::endl;
  }

  int n_iter = 0;
  do {
    ++n_iter;
    if (PRINT_DEBUG_INFO) {
      std::cout << "Iteration #" << n_iter << std::endl;
    }

    // (1a)
    gpuMultiplyKlocalByVec(gpu_data, FEMdata.DIM, FEMdata.elementsCount);
    // (1b)
    float sumElem = gpuDotProduct2(gpu_data.get_p(), gpu_data.get_u(), grid_size);
    // (1c,d)
    float alpha = gamma / sumElem;
    // (2a)
    gpuAddWeighted2(gpu_data.get_x(), gpu_data.get_p(), 1.0f, alpha, grid_size);
    // (2b)
    gpuAddWeighted2(gpu_data.get_r(), gpu_data.get_u(), 1.0f, -1.0f * alpha, grid_size);
    // (3)
    gpuDivideByElementwise(gpu_data.get_r(), gpu_data.get_m(), gpu_data.get_z(), grid_size);
    // (4)
    gpuReductionWithMaskAndTransform2(gpu_data.get_z(), gpu_data.get_mask(), grid_size, gpu_data.get_s(), n_gl_dofs);
    gamma_new = gpuDotProduct2(gpu_data.get_r(), gpu_data.get_s(), grid_size);

    if (PRINT_DEBUG_INFO) {
      // Verbose
      // -----------------------------------------------------------------------
      std::cout << "alpha (gamma / sumElem)\t= " << alpha << std::endl;
      std::cout << "alpha numerator (gamma)\t= " << gamma << std::endl;
      std::cout << "alpha denominator\t= " << sumElem << std::endl;
      // -----------------------------------------------------------------------
    }
    // (5)
    if (gamma_new < eps * gamma0)
      break;

    // (6)
    gpuAddWeighted2(gpu_data.get_p(), gpu_data.get_s(), gamma_new / gamma, 1.0f, grid_size);
    if (PRINT_DEBUG_INFO) {
      // Verbose
      // -----------------------------------------------------------------------
      std::cout << "beta\t\t\t= " << gamma_new / gamma << std::endl;
      std::cout << "gamma_new\t\t= " << gamma_new << std::endl;
      std::cout << std::endl;
      // -----------------------------------------------------------------------
    }
    gamma = gamma_new;
  } while (1);
  if (doAssemblyRes) {
    gpuReductionWithMask2(gpu_data.get_x(), gpu_data.get_mask(), grid_size, gpu_data.get_temp_res());
    gpuDivideByElementwise(gpu_data.get_temp_res(), gpu_data.get_n_adjelem(), gpu_data.get_temp_res(), n_gl_dofs);
    gpuCopyDeviceToHost(gpu_data.get_temp_res(), res.get_data(), n_gl_dofs);
  } else {
    gpuCopyDeviceToHost(gpu_data.get_x(), res.get_data(), grid_size);
  }

}

void gpuCalculateFEM_DYN(FEMdataKeeper &FEMdata, float rho, float damping_alpha, float damping_beta, float endtime, float dt, float beta1, float beta2, bool PRINT_DEBUG_INFO) {
  if (beta2 == 0.0f)    // explicit scheme
    gpuCalculateFEM_explicit(FEMdata, rho, damping_alpha, damping_beta, endtime, dt, beta1, PRINT_DEBUG_INFO);
  else {                // implicit scheme
    gpuCalculateFEM_implicit(FEMdata, rho, damping_alpha, damping_beta, endtime, dt, beta1, beta2, PRINT_DEBUG_INFO);
  }
}

void gpuCalculateFEM_implicit(FEMdataKeeper &FEMdata, float rho, float damping_alpha, float damping_beta, float endtime, float dt, float beta1, float beta2, bool PRINT_DEBUG_INFO) {
  bool is_damping = !((damping_alpha == 0.0f) && (damping_beta == 0.0f));
  if (is_damping)
    gpuCalculateFEM_implicit_DYN_DAMP(FEMdata, rho, damping_alpha, damping_beta, endtime, dt, beta1, beta2, PRINT_DEBUG_INFO);
  else
    gpuCalculateFEM_implicit_DYN(FEMdata, rho, endtime, dt, beta1, beta2, PRINT_DEBUG_INFO);
}

void gpuCalculateFEM_implicit_DYN(FEMdataKeeper &FEMdata, float rho, float endtime, float dt, float beta1, float beta2, bool PRINT_DEBUG_INFO) {
  // Zienkiewicz, Taylor, Zhu "The Finite Element Method: Its Basis and Fundamentals" 6th edition 17.3.3 GN22 (page 608)
  CheckRunTime(__func__)
  //assert(beta2 >= beta1 && beta1 >= 0.5f);

  int n_elems = FEMdata.elementsCount;
  int n_gl_dofs = FEMdata.nodesCount * FEMdata.DIM;
  int grid_size = 3 * FEMdata.DIM * n_elems;
  bool doAssemblyRes = false;
  bool isLumped = false;
  float eps_PCG = 1e-4f;

  for (std::vector<BoundaryEdge>::iterator it = FEMdata.boundary.begin(); it != FEMdata.boundary.end(); ++it) {
    // ERROR:
    // FEMdata.elements[it->adj_elem1].CalculateFlocal(*it, FEMdata.nodes, FEMdata.pressure[num++]);
    if (FEMdata.DIM == 2)
      FEMdata.elements[it->adj_elem1].CalculateFlocal2D(*it, FEMdata.nodes, 0.0f);
    else if (FEMdata.DIM == 3)
      FEMdata.elements[it->adj_elem1].CalculateFlocal3D(*it, FEMdata.nodes, 0.0f);
  }

  std::unordered_map <int, std::vector<int>> nodeAdjElem;
  CalculateNodeAdjElem(FEMdata, nodeAdjElem);

  AssignLoadElement(FEMdata, nodeAdjElem);
//  ApplyLoads_EbE(FEMdata);

  ApplyConstraints_EbE(FEMdata);

//    for (int eIdx = 0; eIdx < n_elems; ++eIdx) {
//        Element *elem = &FEMdata.elements[eIdx];
//        elem->CalculateKlocal(FEMdata.D, FEMdata.nodes);
//        elem->CalculateMlocal(rho, FEMdata.nodes, true); // Lumping on
//        if (is_damping) elem->CalculateClocal(damping_alpha, damping_beta);
//    }

  gpuDataKeeper_DYN gpu_data(FEMdata.DIM, FEMdata.elementsCount, FEMdata.nodesCount, doAssemblyRes, isLumped);
  gpu_data.set_SLAU_matrix_coefs(1.0f, 0.5f * beta2 * dt*dt);
  copyElementsAndFlocals(FEMdata, gpu_data);
  gpuGenerateMask(gpu_data, FEMdata.DIM, FEMdata.elementsCount);
  gpuCountNAdjElem(gpu_data, grid_size);
  gpuCalculateKlocal2(gpu_data, FEMdata);
  gpu_data.copyBmatrixToHost(FEMdata.all_B.get_data());
  gpuCalculateMlocal(gpu_data, FEMdata, rho);
  gpuCalculateDiag(gpu_data, n_elems);

  int endnt = static_cast<int>(endtime / dt);
  int nt = 1;
  while (nt <= endnt) {
      float t = nt*dt;
      if (PRINT_DEBUG_INFO) {
        std::cout << "======= Time iteration #" << nt << "/" << endnt << " =======" << std::endl;
        std::cout << "=========== Time " << t << " ===========" << std::endl << std::endl;
      }

      gpuAddWeighted2(gpu_data.get_displ(), gpu_data.get_vel(), 1.0f, dt, grid_size);
      gpuAddWeighted2(gpu_data.get_displ(), gpu_data.get_x(), 1.0f, 0.5f * (1.0f - beta2) * dt*dt, grid_size);
      gpuAddWeighted2(gpu_data.get_vel(), gpu_data.get_x(), 1.0f, (1.0f - beta1) * dt, grid_size);
      gpu_data.setZeroVec();

      gpuMultiplyMatrixByVec(gpu_data.get_Klocals(), gpu_data.get_displ(), gpu_data.get_r(), FEMdata.elementsCount);

      // update F
      for (int beIdx = 0; beIdx < FEMdata.boundaryEdgesCount; ++beIdx) {
          BoundaryEdge bedge = FEMdata.boundary[beIdx];
          int eIdx = bedge.adj_elem1;
          if (FEMdata.DIM == 2)
            FEMdata.elements[eIdx].CalculateFlocal2D(bedge, FEMdata.nodes, t);
          else if (FEMdata.DIM == 3)
            FEMdata.elements[eIdx].CalculateFlocal3D(bedge, FEMdata.nodes, t);
          // ToDO: ADD APPLY CONSTRAINTS FOR Flocal!
      }
      copyFlocals(FEMdata, gpu_data);

      gpuAddWeighted2(gpu_data.get_r(), gpu_data.get_Flocals(), -1.0f, 1.0f, grid_size);

      // Think how not to use loadVectors! Too dificult!
      std::unordered_map <int, MyArray> loadVectors;
      loadVectors.clear();    // in order to GetMapElement2Loadvector, because if not cleared the values are added
                              // instead of assigned. See if-statement in for-loop in the function's body
      GetMapElement2Loadvector(FEMdata, loadVectors, t);
      copyLoads(gpu_data, loadVectors, FEMdata.DIM, n_elems);

      gpuAdd(gpu_data.get_r(), gpu_data.get_loads(), grid_size);

      gpuPCG_EbE_vec_DYN(FEMdata, gpu_data, doAssemblyRes, eps_PCG, PRINT_DEBUG_INFO);

      gpuAddWeighted2(gpu_data.get_vel(), gpu_data.get_x(), 1.0f, beta1 * dt, grid_size);
      gpuAddWeighted2(gpu_data.get_displ(), gpu_data.get_x(), 1.0f, 0.5f * beta2 * dt*dt, grid_size);


      if (nt % 1 == 0) {
        std::cout << "DEBUG: calling writeSnapshot\n";
        writeSnapshot(t, 21, grid_size, n_gl_dofs, FEMdata, gpu_data);
      }

      ++nt;
      if (PRINT_DEBUG_INFO) {
        std::cout << std::endl;
      }
  }

  gpuReductionWithMask2(gpu_data.get_displ(), gpu_data.get_mask(), grid_size, gpu_data.get_displ_global());
  gpuDivide(gpu_data.get_displ_global(), gpu_data.get_n_adjelem(), n_gl_dofs);        // displ = displ./n_adjelem

  gpuCopyDeviceToHost(gpu_data.get_displ_global(), FEMdata.displacements.get_data(), n_gl_dofs);
}

void gpuCalculateFEM_implicit_DYN_DAMP(FEMdataKeeper &FEMdata, float rho, float damping_alpha, float damping_beta, float endtime, float dt, float beta1, float beta2, bool PRINT_DEBUG_INFO) {
  // Zienkiewicz, Taylor, Zhu "The Finite Element Method: Its Basis and Fundamentals" 6th edition 17.3.3 GN22 (page 608)
  CheckRunTime(__func__)
  //assert(beta2 >= beta1 && beta1 >= 0.5f);

  int n_elems = FEMdata.elementsCount;
  int n_gl_dofs = FEMdata.nodesCount * FEMdata.DIM;
  int grid_size = 3 * FEMdata.DIM * n_elems;
  bool doAssemblyRes = false;
  bool isLumped = false;
  float eps_PCG   = 1e-4f;
  float eps_relax = 1e-4f;
  bool is_relax;

  for (std::vector<BoundaryEdge>::iterator it = FEMdata.boundary.begin(); it != FEMdata.boundary.end(); ++it) {
    // ERROR:
    // FEMdata.elements[it->adj_elem1].CalculateFlocal(*it, FEMdata.nodes, FEMdata.pressure[num++]);
    if (FEMdata.DIM == 2)
      FEMdata.elements[it->adj_elem1].CalculateFlocal2D(*it, FEMdata.nodes, 0.0f);
    else if (FEMdata.DIM == 3)
      FEMdata.elements[it->adj_elem1].CalculateFlocal3D(*it, FEMdata.nodes, 0.0f);
  }

  std::unordered_map <int, std::vector<int>> nodeAdjElem;
  CalculateNodeAdjElem(FEMdata, nodeAdjElem);

  AssignLoadElement(FEMdata, nodeAdjElem);
//  ApplyLoads_EbE(FEMdata);

  ApplyConstraints_EbE(FEMdata);

//    for (int eIdx = 0; eIdx < n_elems; ++eIdx) {
//        Element *elem = &FEMdata.elements[eIdx];
//        elem->CalculateKlocal(FEMdata.D, FEMdata.nodes);
//        elem->CalculateMlocal(rho, FEMdata.nodes, true); // Lumping on
//        if (is_damping) elem->CalculateClocal(damping_alpha, damping_beta);
//    }

  gpuDataKeeper_DYN_DAMP gpu_data(FEMdata.DIM, FEMdata.elementsCount, FEMdata.nodesCount, doAssemblyRes, isLumped, damping_alpha, damping_beta);
  gpu_data.set_SLAU_matrix_coefs(1.0f, 0.5f * beta2 * dt*dt, beta1 * dt);
  copyElementsAndFlocals(FEMdata, gpu_data);
  gpuGenerateMask(gpu_data, FEMdata.DIM, FEMdata.elementsCount);
  gpuCountNAdjElem(gpu_data, grid_size);
  gpuCalculateKlocal2(gpu_data, FEMdata);
  gpu_data.copyBmatrixToHost(FEMdata.all_B.get_data());
  gpuCalculateMlocal(gpu_data, FEMdata, rho);
  gpuCalculateDiag_DAMP(gpu_data, n_elems);

  int endnt;
  is_relax = (endtime < 0.0f);
  if (!is_relax) endnt = static_cast<int>(endtime / dt);
  int nt = 1;
  float cnorm_acc, cnorm_vel;
  do {
      float t = nt*dt;
      if (PRINT_DEBUG_INFO) {
        std::cout << "======= Time iteration #" << nt;
        if (!is_relax) cout << "/" << endnt << " =======";
        std::cout << " =======";
        std::cout << "\n=========== Time " << t << " ===========\n\n";
      }

      gpuAddWeighted2(gpu_data.get_displ(), gpu_data.get_vel(), 1.0f, dt, grid_size);
      gpuAddWeighted2(gpu_data.get_displ(), gpu_data.get_x(), 1.0f, 0.5f * (1.0f - beta2) *  dt*dt, grid_size);
      gpuAddWeighted2(gpu_data.get_vel(), gpu_data.get_x(), 1.0f, (1.0f - beta1) * dt, grid_size);
      gpu_data.setZeroVec();

      gpuMultiplyMatrixByVec(gpu_data.get_Klocals(), gpu_data.get_displ(), gpu_data.get_r(), FEMdata.elementsCount);
      //if (is_damping) {
      gpuMultiplyClocalByVec(gpu_data, gpu_data.get_vel(), gpu_data.get_tmp(), FEMdata.elementsCount);
      gpuAdd(gpu_data.get_r(), gpu_data.get_tmp(), grid_size);
      //}

      // update F
      for (int beIdx = 0; beIdx < FEMdata.boundaryEdgesCount; ++beIdx) {
          BoundaryEdge bedge = FEMdata.boundary[beIdx];
          int eIdx = bedge.adj_elem1;
          if (FEMdata.DIM == 2)
            FEMdata.elements[eIdx].CalculateFlocal2D(bedge, FEMdata.nodes, t);
          else if (FEMdata.DIM == 3)
            FEMdata.elements[eIdx].CalculateFlocal3D(bedge, FEMdata.nodes, t);
          // ToDO: ADD APPLY CONSTRAINTS FOR Flocal!
      }
      copyFlocals(FEMdata, gpu_data);

      gpuAddWeighted2(gpu_data.get_r(), gpu_data.get_Flocals(), -1.0f, 1.0f, grid_size);

      // Think how not to use loadVectors! Too dificult!
      std::unordered_map <int, MyArray> loadVectors;
      loadVectors.clear();    // in order to GetMapElement2Loadvector, because if not cleared the values are added
                              // instead of assigned. See if-statement in for-loop in the function's body
      GetMapElement2Loadvector(FEMdata, loadVectors, t);
      copyLoads(gpu_data, loadVectors, FEMdata.DIM, n_elems);

      gpuAdd(gpu_data.get_r(), gpu_data.get_loads(), grid_size);

      gpuPCG_EbE_vec_DYN_DAMP(FEMdata, gpu_data, doAssemblyRes, eps_PCG, PRINT_DEBUG_INFO);

      gpuAddWeighted2(gpu_data.get_vel(), gpu_data.get_x(), 1.0f, beta1 * dt, grid_size);
      gpuAddWeighted2(gpu_data.get_displ(), gpu_data.get_x(), 1.0f, 0.5f * beta2 * dt*dt, grid_size);


      ++nt;

      if (is_relax) {
        cnorm_acc = gpuCNorm(gpu_data.get_x(), grid_size);
        cnorm_vel = gpuCNorm(gpu_data.get_vel(), grid_size);
        if (PRINT_DEBUG_INFO) {
          std::cout << "C-norm acc = " << cnorm_acc << "\nC-norm vel = " << cnorm_vel << std::endl;
          std::cout << std::endl;
        }
        if ((cnorm_vel < eps_relax)) break; // && (cnorm_acc < eps_relax)
      } else {
        if (PRINT_DEBUG_INFO) {
          std::cout << std::endl;
        }

        if (nt > endnt) break;
      }
  } while (true);

  gpuReductionWithMask2(gpu_data.get_displ(), gpu_data.get_mask(), grid_size, gpu_data.get_displ_global());
  gpuDivide(gpu_data.get_displ_global(), gpu_data.get_n_adjelem(), n_gl_dofs);        // displ = displ./n_adjelem

  gpuCopyDeviceToHost(gpu_data.get_displ_global(), FEMdata.displacements.get_data(), n_gl_dofs);
}

void gpuCalculateFEM_explicit(FEMdataKeeper &FEMdata, float rho, float damping_alpha, float damping_beta, float endtime, float dt, float beta1, bool PRINT_DEBUG_INFO) {
    bool is_damping = !((damping_alpha == 0.0f) && (damping_beta == 0.0f));
    if (is_damping)
      gpuCalculateFEM_explicit_DYN_DAMP(FEMdata, rho, damping_alpha, damping_beta, endtime, dt, beta1, PRINT_DEBUG_INFO);
    else
      gpuCalculateFEM_explicit_DYN(FEMdata, rho, endtime, dt, beta1, PRINT_DEBUG_INFO);
}

void gpuCalculateFEM_explicit_DYN(FEMdataKeeper &FEMdata, float rho, float endtime, float dt, float beta1, bool PRINT_DEBUG_INFO) {
  // Zienkiewicz, Taylor, Zhu "The Finite Element Method: Its Basis and Fundamentals" 6th edition 17.3.3 GN22 (page 608)
  CheckRunTime(__func__)

  int n_elems = FEMdata.elementsCount;
  int n_gl_dofs = FEMdata.nodesCount * FEMdata.DIM;
  int grid_size = 3 * FEMdata.DIM * n_elems;
  bool doAssemblyRes = false;
  bool isLumped = true;
  float eps_relax = 1e-4f;
  bool is_relax;

  for (std::vector<BoundaryEdge>::iterator it = FEMdata.boundary.begin(); it != FEMdata.boundary.end(); ++it) {
    // ERROR:
    // FEMdata.elements[it->adj_elem1].CalculateFlocal(*it, FEMdata.nodes, FEMdata.pressure[num++]);
    if (FEMdata.DIM == 2)
      FEMdata.elements[it->adj_elem1].CalculateFlocal2D(*it, FEMdata.nodes, 0.0f);
    else if (FEMdata.DIM == 3)
      FEMdata.elements[it->adj_elem1].CalculateFlocal3D(*it, FEMdata.nodes, 0.0f);
  }

  std::unordered_map <int, std::vector<int>> nodeAdjElem;
  CalculateNodeAdjElem(FEMdata, nodeAdjElem);

  AssignLoadElement(FEMdata, nodeAdjElem);
//  ApplyLoads_EbE(FEMdata);

  ApplyConstraints_EbE(FEMdata);

//    for (int eIdx = 0; eIdx < n_elems; ++eIdx) {
//        Element *elem = &FEMdata.elements[eIdx];
//        elem->CalculateKlocal(FEMdata.D, FEMdata.nodes);
//        elem->CalculateMlocal(rho, FEMdata.nodes, true); // Lumping on
//        if (is_damping) elem->CalculateClocal(damping_alpha, damping_beta);
//    }

  gpuDataKeeper_DYN gpu_data(FEMdata.DIM, FEMdata.elementsCount, FEMdata.nodesCount, doAssemblyRes, isLumped);
  copyElementsAndFlocals(FEMdata, gpu_data);
  gpuGenerateMask(gpu_data, FEMdata.DIM, FEMdata.elementsCount);
  gpuCountNAdjElem(gpu_data, grid_size);
  gpuCalculateKlocal2(gpu_data, FEMdata);
  gpu_data.copyBmatrixToHost(FEMdata.all_B.get_data());
  gpuCalculateMlocal(gpu_data, FEMdata, rho);

  int endnt;
  is_relax = (endtime < 0.0f);
  if (!is_relax) endnt = static_cast<int>(endtime / dt);
  int nt = 1;
  float cnorm_acc, cnorm_vel;
  do {
    float t = nt*dt;
    if (PRINT_DEBUG_INFO) {
      std::cout << "======= Time iteration #" << nt;
      if (!is_relax) cout << "/" << endnt;
      std::cout << " =======";
      std::cout << "\n=========== Time " << t << " ===========\n\n";
    }

      gpuAddWeighted2(gpu_data.get_displ(), gpu_data.get_vel(), 1.0f, dt, grid_size);
      gpuAddWeighted2(gpu_data.get_displ(), gpu_data.get_x(), 1.0f, 0.5f * dt*dt, grid_size);
      gpuAddWeighted2(gpu_data.get_vel(), gpu_data.get_x(), 1.0f, (1.0f - beta1) * dt, grid_size);

      gpuMultiplyMatrixByVec(gpu_data.get_Klocals(), gpu_data.get_displ(), gpu_data.get_r(), FEMdata.elementsCount);

      // update F
      for (int beIdx = 0; beIdx < FEMdata.boundaryEdgesCount; ++beIdx) {
          BoundaryEdge bedge = FEMdata.boundary[beIdx];
          int eIdx = bedge.adj_elem1;
          if (FEMdata.DIM == 2)
            FEMdata.elements[eIdx].CalculateFlocal2D(bedge, FEMdata.nodes, t);
          else if (FEMdata.DIM == 3)
            FEMdata.elements[eIdx].CalculateFlocal3D(bedge, FEMdata.nodes, t);
          // ToDO: ADD APPLY CONSTRAINTS FOR Flocal!
      }
      copyFlocals(FEMdata, gpu_data);

      gpuAddWeighted2(gpu_data.get_r(), gpu_data.get_Flocals(), -1.0f, 1.0f, grid_size);

      // Think how not to use loadVectors! Too dificult!
      std::unordered_map <int, MyArray> loadVectors;
      loadVectors.clear();    // in order to GetMapElement2Loadvector, because if not cleared the values are added
                              // instead of assigned. See if-statement in for-loop in the function's body
      GetMapElement2Loadvector(FEMdata, loadVectors, t);
      copyLoads(gpu_data, loadVectors, FEMdata.DIM, n_elems);

      gpuAdd(gpu_data.get_r(), gpu_data.get_loads(), grid_size);

      gpuSolveDiag(gpu_data.get_diagM(), gpu_data.get_r(),
                   gpu_data.get_x(), gpu_data.get_mask(),
                   gpu_data.get_n_adjelem(), grid_size, n_gl_dofs, doAssemblyRes);

      gpuAddWeighted2(gpu_data.get_vel(), gpu_data.get_x(), 1.0f, beta1 * dt, grid_size);

      if (nt % 10) {
        writeSnapshot(t, 20, grid_size, n_gl_dofs, FEMdata, gpu_data);
      }

      ++nt;

      if (is_relax) {
        cnorm_acc = gpuCNorm(gpu_data.get_x(), grid_size);
        cnorm_vel = gpuCNorm(gpu_data.get_vel(), grid_size);
        if (PRINT_DEBUG_INFO) {
          std::cout << "C-norm acc = " << cnorm_acc << "\nC-norm vel = " << cnorm_vel << std::endl;
          std::cout << std::endl;
        }
        if ((cnorm_vel < eps_relax)) break; // && (cnorm_acc < eps_relax)
      } else {
        if (PRINT_DEBUG_INFO) {
          std::cout << std::endl;
        }

        if (nt > endnt) break;
      }
  } while (true);

  gpuReductionWithMask2(gpu_data.get_displ(), gpu_data.get_mask(), grid_size, gpu_data.get_displ_global());
  gpuDivide(gpu_data.get_displ_global(), gpu_data.get_n_adjelem(), n_gl_dofs);        // displ = displ./n_adjelem

  gpuCopyDeviceToHost(gpu_data.get_displ_global(), FEMdata.displacements.get_data(), n_gl_dofs);
}

void gpuCalculateFEM_explicit_DYN_DAMP(FEMdataKeeper &FEMdata, float rho, float damping_alpha, float damping_beta, float endtime, float dt, float beta1, bool PRINT_DEBUG_INFO) {
  // Zienkiewicz, Taylor, Zhu "The Finite Element Method: Its Basis and Fundamentals" 6th edition 17.3.3 GN22 (page 608)
  CheckRunTime(__func__)
  assert(damping_beta == 0.0f); // ToDO: rewrite considering this!

  fstream out1;
  out1.open("C:/Users/mokin/Desktop/fem_stuff_test/fem_results/graph/" + FEMdata.name + ".txt", fstream::out);
  out1 << FEMdata.nodesCount << " " << FEMdata.elementsCount << "\n";

  int n_elems = FEMdata.elementsCount;
  int n_gl_dofs = FEMdata.nodesCount * FEMdata.DIM;
  int grid_size = 3 * FEMdata.DIM * n_elems;
  bool doAssemblyRes = false;
  bool isLumped = true;
  float eps_relax;
  bool is_relax;

  for (std::vector<BoundaryEdge>::iterator it = FEMdata.boundary.begin(); it != FEMdata.boundary.end(); ++it) {
    // ERROR:
    // FEMdata.elements[it->adj_elem1].CalculateFlocal(*it, FEMdata.nodesX, FEMdata.nodesY, FEMdata.pressure[num++]);
    if (FEMdata.DIM == 2)
      FEMdata.elements[it->adj_elem1].CalculateFlocal2D(*it, FEMdata.nodes, 0.0f);
    else if (FEMdata.DIM == 3)
      FEMdata.elements[it->adj_elem1].CalculateFlocal3D(*it, FEMdata.nodes, 0.0f);
  }

  std::unordered_map <int, std::vector<int>> nodeAdjElem;
  CalculateNodeAdjElem(FEMdata, nodeAdjElem);

  AssignLoadElement(FEMdata, nodeAdjElem);
//  ApplyLoads_EbE(FEMdata);

  ApplyConstraints_EbE(FEMdata);

//    for (int eIdx = 0; eIdx < n_elems; ++eIdx) {
//        Element *elem = &FEMdata.elements[eIdx];
//        elem->CalculateKlocal(FEMdata.D, FEMdata.nodes);
//        elem->CalculateMlocal(rho, FEMdata.nodes, true); // Lumping on
//        if (is_damping) elem->CalculateClocal(damping_alpha, damping_beta);
//    }

  gpuDataKeeper_DYN_DAMP gpu_data(FEMdata.DIM, FEMdata.elementsCount, FEMdata.nodesCount, doAssemblyRes, isLumped, damping_alpha, damping_beta);
  gpu_data.set_SLAU_matrix_coefs(1.0f, 0.0f, beta1 * dt);
  copyElementsAndFlocals(FEMdata, gpu_data);
  gpuGenerateMask(gpu_data, FEMdata.DIM, FEMdata.elementsCount);
  gpuCountNAdjElem(gpu_data, grid_size);
  gpuCalculateKlocal2(gpu_data, FEMdata);
  gpu_data.copyBmatrixToHost(FEMdata.all_B.get_data());
  gpuCalculateMlocal(gpu_data, FEMdata, rho);
  gpuCalculateDiag_DAMP(gpu_data, n_elems);

  int endnt;
  is_relax = (endtime < 0.0f);
  if (!is_relax) endnt = static_cast<int>(endtime / dt);
  else {
    eps_relax = std::fabs(endtime);
    if (PRINT_DEBUG_INFO)
      std::cout << "eps_relax = " << eps_relax << "\n";

  }
  int nt = 1;
  float cnorm_acc, cnorm_vel;
  do {
      float t = nt*dt;
      if (PRINT_DEBUG_INFO) {
        std::cout << "======= Time iteration #" << nt;
        if (!is_relax) cout << "/" << endnt;
        std::cout << " =======";
        std::cout << "\n=========== Time " << t << " ===========\n\n";
      }

      gpuAddWeighted2(gpu_data.get_displ(), gpu_data.get_vel(), 1.0f, dt, grid_size);
      gpuAddWeighted2(gpu_data.get_displ(), gpu_data.get_x(), 1.0f, 0.5f * dt*dt, grid_size);
      gpuAddWeighted2(gpu_data.get_vel(), gpu_data.get_x(), 1.0f, (1.0f - beta1) * dt, grid_size);
      //gpu_data.setZeroVec();

      gpuMultiplyMatrixByVec(gpu_data.get_Klocals(), gpu_data.get_displ(), gpu_data.get_r(), FEMdata.elementsCount);
      //if (is_damping) {
      gpuMultiplyClocalByVec(gpu_data, gpu_data.get_vel(), gpu_data.get_tmp(), FEMdata.elementsCount);
      gpuAdd(gpu_data.get_r(), gpu_data.get_tmp(), grid_size);
      //}

      // update F
      for (int beIdx = 0; beIdx < FEMdata.boundaryEdgesCount; ++beIdx) {
          BoundaryEdge bedge = FEMdata.boundary[beIdx];
          int eIdx = bedge.adj_elem1;
          if (FEMdata.DIM == 2)
            FEMdata.elements[eIdx].CalculateFlocal2D(bedge, FEMdata.nodes, t);
          else if (FEMdata.DIM == 3)
            FEMdata.elements[eIdx].CalculateFlocal3D(bedge, FEMdata.nodes, t);
          // ToDO: ADD APPLY CONSTRAINTS FOR Flocal!
      }
      copyFlocals(FEMdata, gpu_data);

      gpuAddWeighted2(gpu_data.get_r(), gpu_data.get_Flocals(), -1.0f, 1.0f, grid_size);

      // Think how not to use loadVectors! Too dificult!
      std::unordered_map <int, MyArray> loadVectors;
      loadVectors.clear();    // in order to GetMapElement2Loadvector, because if not cleared the values are added
                              // instead of assigned. See if-statement in for-loop in the function's body
      GetMapElement2Loadvector(FEMdata, loadVectors, t);
      copyLoads(gpu_data, loadVectors, FEMdata.DIM, n_elems);

      gpuAdd(gpu_data.get_r(), gpu_data.get_loads(), grid_size);

      gpuSolveDiag(gpu_data.get_diag(), gpu_data.get_r(),
                   gpu_data.get_x(), gpu_data.get_mask(),
                   gpu_data.get_n_adjelem(), grid_size, n_gl_dofs, doAssemblyRes);

      gpuAddWeighted2(gpu_data.get_vel(), gpu_data.get_x(), 1.0f, beta1 * dt, grid_size);


      ++nt;

      if (is_relax) {
        cnorm_acc = gpuCNorm(gpu_data.get_x(), grid_size);
        cnorm_vel = gpuCNorm(gpu_data.get_vel(), grid_size);

        out1 << cnorm_vel << " " << nt * dt << "\n";

        if (PRINT_DEBUG_INFO) {
          std::cout << "C-norm acc = " << cnorm_acc << "\nC-norm vel = " << cnorm_vel << std::endl;
          std::cout << std::endl;
        }
        if ((cnorm_vel < eps_relax)) break; // && (cnorm_acc < eps_relax)
      } else {
        if (PRINT_DEBUG_INFO) {
          std::cout << std::endl;
        }

        if (nt > endnt) break;
      }
  } while (true);

  gpuReductionWithMask2(gpu_data.get_displ(), gpu_data.get_mask(), grid_size, gpu_data.get_displ_global());
  gpuDivide(gpu_data.get_displ_global(), gpu_data.get_n_adjelem(), n_gl_dofs);        // displ = displ./n_adjelem

  gpuCopyDeviceToHost(gpu_data.get_displ_global(), FEMdata.displacements.get_data(), n_gl_dofs);
}

void gpuPCG_EbE_vec_DYN(FEMdataKeeper &FEMdata, gpuDataKeeper_DYN &gpu_data, bool doAssemblyRes,  float eps, bool PRINT_DEBUG_INFO) {
  CheckRunTime(__func__)
  int n_elems  = FEMdata.elementsCount;
  int n_gl_dofs = FEMdata.nodesCount * FEMdata.DIM;
  int grid_size = 3 * FEMdata.DIM * n_elems;

  // (0a)
  // Initialize vector r^(e)
  //gpuCopyDeviceToDevice(gpu_data.get_b(), gpu_data.get_r(), grid_size);
  // (0b)
  gpuReductionWithMaskAndTransform2(gpu_data.get_diag(), gpu_data.get_mask(), grid_size, gpu_data.get_m(), n_gl_dofs);
  // (0c)
  gpuDivideByElementwise(gpu_data.get_r(), gpu_data.get_m(), gpu_data.get_z(), grid_size);
  // (0d)
  gpuReductionWithMaskAndTransform2(gpu_data.get_z(), gpu_data.get_mask(), grid_size, gpu_data.get_s(), n_gl_dofs);
  float gamma0 = gpuDotProduct2(gpu_data.get_r(), gpu_data.get_s(), grid_size);
  float gamma = gamma0;
  float gamma_new = 0.0f;

  // (0e)
  gpuCopyDeviceToDevice(gpu_data.get_s(), gpu_data.get_p(), grid_size);

  if (PRINT_DEBUG_INFO) {
    //std::cout.precision(16);
    std::cout << "gamma0\t\t\t= " << gamma0 << std::endl << std::endl;
  }

  int n_iter = 0;
  do {
    ++n_iter;
    if (PRINT_DEBUG_INFO) {
      std::cout << "Iteration #" << n_iter << std::endl;
    }

    // (1a)
    //gpuMultiplyKlocalByVec(gpu_data, FEMdata.elementsCount);
    gpuMultiplyAlocalByVec(gpu_data, FEMdata.elementsCount);
    // (1b)
    float sumElem = gpuDotProduct2(gpu_data.get_p(), gpu_data.get_u(), grid_size);
    // (1c,d)
    float alpha = gamma / sumElem;
    // (2a)
    gpuAddWeighted2(gpu_data.get_x(), gpu_data.get_p(), 1.0f, alpha, grid_size);
    // (2b)
    gpuAddWeighted2(gpu_data.get_r(), gpu_data.get_u(), 1.0f, -1.0f * alpha, grid_size);
    // (3)
    gpuDivideByElementwise(gpu_data.get_r(), gpu_data.get_m(), gpu_data.get_z(), grid_size);
    // (4)
    gpuReductionWithMaskAndTransform2(gpu_data.get_z(), gpu_data.get_mask(), grid_size, gpu_data.get_s(), n_gl_dofs);
    gamma_new = gpuDotProduct2(gpu_data.get_r(), gpu_data.get_s(), grid_size);

    if (PRINT_DEBUG_INFO) {
      // Verbose
      // -----------------------------------------------------------------------
      std::cout << "alpha (gamma / sumElem)\t= " << alpha << std::endl;
      std::cout << "alpha numerator (gamma)\t= " << gamma << std::endl;
      std::cout << "alpha denominator\t= " << sumElem << std::endl;
      // -----------------------------------------------------------------------
    }
    // (5)
    if (gamma_new < eps * gamma0)
      break;

    // (6)
    gpuAddWeighted2(gpu_data.get_p(), gpu_data.get_s(), gamma_new / gamma, 1.0f, grid_size);
    if (PRINT_DEBUG_INFO) {
      // Verbose
      // -----------------------------------------------------------------------
      std::cout << "beta\t\t\t= " << gamma_new / gamma << std::endl;
      std::cout << "gamma_new\t\t= " << gamma_new << std::endl;
      std::cout << std::endl;
      // -----------------------------------------------------------------------
    }
    gamma = gamma_new;
  } while (1);
  if (doAssemblyRes) {
    gpuReductionWithMask2(gpu_data.get_x(), gpu_data.get_mask(), grid_size, gpu_data.get_temp_res());
    gpuDivideByElementwise(gpu_data.get_temp_res(), gpu_data.get_n_adjelem(), gpu_data.get_temp_res(), n_gl_dofs);
  }

}

void gpuPCG_EbE_vec_DYN_DAMP(FEMdataKeeper &FEMdata, gpuDataKeeper_DYN_DAMP &gpu_data, bool doAssemblyRes,  float eps, bool PRINT_DEBUG_INFO) {
  CheckRunTime(__func__)
  int n_elems  = FEMdata.elementsCount;
  int n_gl_dofs = FEMdata.nodesCount * FEMdata.DIM;
  int grid_size = 3 * FEMdata.DIM * n_elems;

  // (0b)
  gpuReductionWithMaskAndTransform2(gpu_data.get_diag(), gpu_data.get_mask(), grid_size, gpu_data.get_m(), n_gl_dofs);
  // (0c)
  gpuDivideByElementwise(gpu_data.get_r(), gpu_data.get_m(), gpu_data.get_z(), grid_size);
  // (0d)
  gpuReductionWithMaskAndTransform2(gpu_data.get_z(), gpu_data.get_mask(), grid_size, gpu_data.get_s(), n_gl_dofs);
  float gamma0 = gpuDotProduct2(gpu_data.get_r(), gpu_data.get_s(), grid_size);
  float gamma = gamma0;
  float gamma_new = 0.0f;

  // (0e)
  gpuCopyDeviceToDevice(gpu_data.get_s(), gpu_data.get_p(), grid_size);

  if (PRINT_DEBUG_INFO) {
    //std::cout.precision(16);
    std::cout << "gamma0\t\t\t= " << gamma0 << std::endl << std::endl;
  }

  int n_iter = 0;
  do {
    ++n_iter;
    if (PRINT_DEBUG_INFO) {
      std::cout << "Iteration #" << n_iter << std::endl;
    }

    // (1a)
    //gpuMultiplyKlocalByVec(gpu_data, FEMdata.elementsCount);
    gpuMultiplyAlocalByVec_DAMP(gpu_data, FEMdata.elementsCount);
    // (1b)
    float sumElem = gpuDotProduct2(gpu_data.get_p(), gpu_data.get_u(), grid_size);
    // (1c,d)
    float alpha = gamma / sumElem;
    // (2a)
    gpuAddWeighted2(gpu_data.get_x(), gpu_data.get_p(), 1.0f, alpha, grid_size);
    // (2b)
    gpuAddWeighted2(gpu_data.get_r(), gpu_data.get_u(), 1.0f, -1.0f * alpha, grid_size);
    // (3)
    gpuDivideByElementwise(gpu_data.get_r(), gpu_data.get_m(), gpu_data.get_z(), grid_size);
    // (4)
    gpuReductionWithMaskAndTransform2(gpu_data.get_z(), gpu_data.get_mask(), grid_size, gpu_data.get_s(), n_gl_dofs);
    gamma_new = gpuDotProduct2(gpu_data.get_r(), gpu_data.get_s(), grid_size);

    if (PRINT_DEBUG_INFO) {
      // Verbose
      // -----------------------------------------------------------------------
      std::cout << "alpha (gamma / sumElem)\t= " << alpha << std::endl;
      std::cout << "alpha numerator (gamma)\t= " << gamma << std::endl;
      std::cout << "alpha denominator\t= " << sumElem << std::endl;
      // -----------------------------------------------------------------------
    }
    // (5)
    if (gamma_new < eps * gamma0)
      break;

    // (6)
    gpuAddWeighted2(gpu_data.get_p(), gpu_data.get_s(), gamma_new / gamma, 1.0f, grid_size);
    if (PRINT_DEBUG_INFO) {
      // Verbose
      // -----------------------------------------------------------------------
      std::cout << "beta\t\t\t= " << gamma_new / gamma << std::endl;
      std::cout << "gamma_new\t\t= " << gamma_new << std::endl;
      std::cout << std::endl;
      // -----------------------------------------------------------------------
    }
    gamma = gamma_new;
  } while (1);
  if (doAssemblyRes) {
    gpuReductionWithMask2(gpu_data.get_x(), gpu_data.get_mask(), grid_size, gpu_data.get_temp_res());
    gpuDivideByElementwise(gpu_data.get_temp_res(), gpu_data.get_n_adjelem(), gpu_data.get_temp_res(), n_gl_dofs);
  }

}


void CalculateFEM_dyn(FEMdataKeeper &FEMdata, float rho, float damping_alpha, float damping_beta, float endtime, float dt, bool PRINT_DEBUG_INFO) {
  // Zienkiewicz, Taylor, Zhu "The Finite Element Method: Its Basis and Fundamentals" 6th edition 17.3.3 GN22 (page 608)
  CheckRunTime(__func__)
  const float beta1 = 0.55f;
  const float beta2 = 0.5f;
  // beta2 = 0.0 -- explicit scheme (assuming both M and C are diagonal)
  // implicit scheme: beta1 >= beta2 >= 1/2
  assert(beta2 == 0.0f || beta2 >= beta1);

  for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
    it->CalculateKlocal(FEMdata.D, FEMdata.nodes);
    it->CalculateMlocal(rho, FEMdata.nodes, false);
    it->CalculateClocal(damping_alpha, damping_beta);
  }

  for (std::vector<BoundaryEdge>::iterator it = FEMdata.boundary.begin(); it != FEMdata.boundary.end(); ++it) {
    if (FEMdata.DIM == 2)
      FEMdata.elements[it->adj_elem1].CalculateFlocal2D(*it, FEMdata.nodes, 0.0f);
    else if (FEMdata.DIM == 3)
      FEMdata.elements[it->adj_elem1].CalculateFlocal3D(*it, FEMdata.nodes, 0.0f);
  }

  std::unordered_map <int, std::vector<int>> nodeAdjElem;
  CalculateNodeAdjElem(FEMdata, nodeAdjElem);

  ApplyConstraints_EbE(FEMdata);

  for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
    it->Alocal = it->Mlocal.weightedSum(it->Clocal, 1.0f, beta1 * dt);
    it->Alocal = it->Alocal.weightedSum(it->Klocal, 1.0f, 0.5f * beta2 * dt*dt);
  }

  AssignLoadElement(FEMdata, nodeAdjElem);

  std::unordered_map <int, MyArray> loadVectors;
  //GetMapElement2Loadvector(FEMdata, loadVectors_old, 0.0f);

  // Set initial condition. Assumed that zero. ToDO: add initial conditions to prepared_meshes

  int endnt = static_cast<int>(endtime / dt);
  int nt = 1;
  while (nt <= endnt) {
    float t = nt*dt;
    std::cout << "======= Time iteration #" << nt << "/" << endnt << " =======" << std::endl;
    std::cout << "=========== Time " << t << " ===========" << std::endl << std::endl;
    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
      it->x_pred = it->x;
      it->x_pred.add_weighted(it->vel, 1.0, dt);
      it->x_pred.add_weighted(it->res, 1.0f, 0.5f * (1 - beta2) * dt*dt);

      it->vel_pred = it->vel;
      it->vel_pred.add_weighted(it->res, 1.0f, (1 - beta1) * dt);

      it->blocal = it->Clocal.Product(it->vel_pred);
      it->blocal.add(it->Klocal.Product(it->x_pred));
      it->blocal.scale(-1.0f);
      it->blocal.add(it->Flocal);  // Seems to be the same (that is, Flocal), see formula (17.15)
    }

    GetMapElement2Loadvector(FEMdata, loadVectors, t);
    for (auto& it: loadVectors) {
      FEMdata.elements[it.first].blocal.add(loadVectors[it.first]);
    }

    PCG_EbE(FEMdata, nodeAdjElem, 1e-10f, PRINT_DEBUG_INFO);
    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
      it->x = it->x_pred;
      it->x.add_weighted(it->res, 1.0f, 0.5f * beta2 * dt*dt);
      it->vel = it->vel_pred;
      it->vel.add_weighted(it->res, 1.0f, beta1 * dt);
    }

    ++nt;
    std::cout << std::endl;
  }

  AssemblyX(FEMdata, nodeAdjElem);
}

void CalculateFEM_dyn_vec(FEMdataKeeper &FEMdata, float rho, float damping_alpha, float damping_beta, float endtime, float dt, float beta1, float beta2, bool PRINT_DEBUG_INFO) {
  // Zienkiewicz, Taylor, Zhu "The Finite Element Method: Its Basis and Fundamentals" 6th edition 17.3.3 GN22 (page 608)
  CheckRunTime(__func__)
  //const float beta1 = 0.55f; //0.55f;
  //const float beta2 = 0.5f; //0.5f;
  // beta2 = 0.0 -- explicit scheme (assuming both M and C are diagonal -- make sure to lump mass matrix!)
  // implicit scheme: beta1 >= beta2 >= 1/2
  int DIM = FEMdata.DIM;
  assert(beta2 == 0.0f || (beta1 >= beta2 && beta2 >= 0.5f));
  bool isLumping = (beta2 == 0.0f);

  for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
    it->CalculateKlocal(FEMdata.D, FEMdata.nodes);
    it->CalculateMlocal(rho, FEMdata.nodes, isLumping);
    it->CalculateClocal(damping_alpha, damping_beta);
  }

  for (std::vector<BoundaryEdge>::iterator it = FEMdata.boundary.begin(); it != FEMdata.boundary.end(); ++it) {
    if (DIM == 2)
      FEMdata.elements[it->adj_elem1].CalculateFlocal2D(*it, FEMdata.nodes, 0.0f);
    else if (DIM == 3)
      FEMdata.elements[it->adj_elem1].CalculateFlocal3D(*it, FEMdata.nodes, 0.0f);
  }

  std::unordered_map <int, std::vector<int>> nodeAdjElem;
  CalculateNodeAdjElem(FEMdata, nodeAdjElem);

  ApplyConstraints_EbE(FEMdata);

  for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
    it->Alocal = it->Mlocal.weightedSum(it->Clocal, 1.0f, beta1 * dt);
    it->Alocal = it->Mlocal;
    it->Alocal = it->Alocal.weightedSum(it->Klocal, 1.0f, 0.5f * beta2 * dt*dt);
  }

  AssignLoadElement(FEMdata, nodeAdjElem);

  //std::unordered_map <int, MyArray> loadVectors;

  // Set initial condition. Assumed that zero. ToDO: add initial conditions to prepared_meshes

  // ----------------------------------------------------------------------------------------------------------------
  int n_elems  = FEMdata.elementsCount;
  int n_gl_dofs = FEMdata.nodesCount * DIM;

  int grid_size = 3 * DIM * n_elems;
  MyArray x(grid_size), x_pred(grid_size), vel(grid_size),
      vel_pred(grid_size), b(grid_size), res(grid_size),
      F(grid_size);
  MyArray n_adjelem(n_gl_dofs);
  SparseMatrixCOO Q (grid_size);

  int nonzeroK = 0, nonzeroA = 0, nonzeroC = 0;
  for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
    nonzeroK += it->Klocal.CountNonzero();
    nonzeroA += it->Alocal.CountNonzero();
    nonzeroC += it->Clocal.CountNonzero();
  }

  FEMdata.nonzeroMatrixNumCount = nonzeroA;

  SparseMatrixCOO K(nonzeroK), C(nonzeroC);
  for (int eIdx = 0; eIdx < n_elems; ++eIdx) {
    Element elem = FEMdata.elements[eIdx];

    Q.write_value(3 * DIM * eIdx + 0 * DIM + 0, elem.nodesIds[0] * DIM + 0, 1);
    Q.write_value(3 * DIM * eIdx + 0 * DIM + 1, elem.nodesIds[0] * DIM + 1, 1);
    Q.write_value(3 * DIM * eIdx + 1 * DIM + 0, elem.nodesIds[1] * DIM + 0, 1);
    Q.write_value(3 * DIM * eIdx + 1 * DIM + 1, elem.nodesIds[1] * DIM + 1, 1);
    Q.write_value(3 * DIM * eIdx + 2 * DIM + 0, elem.nodesIds[2] * DIM + 0, 1);
    Q.write_value(3 * DIM * eIdx + 2 * DIM + 1, elem.nodesIds[2] * DIM + 1, 1);

    n_adjelem[elem.nodesIds[0] * DIM + 0] += 1.0f;
    n_adjelem[elem.nodesIds[0] * DIM + 1] += 1.0f;
    n_adjelem[elem.nodesIds[1] * DIM + 0] += 1.0f;
    n_adjelem[elem.nodesIds[1] * DIM + 1] += 1.0f;
    n_adjelem[elem.nodesIds[2] * DIM + 0] += 1.0f;
    n_adjelem[elem.nodesIds[2] * DIM + 1] += 1.0f;

    for (int xlocal = 0; xlocal < 3 * DIM; ++xlocal) {
      for (int ylocal = 0; ylocal < 3 * DIM; ++ylocal) {
        K.write_value(3 * DIM * eIdx + xlocal, 3 * DIM * eIdx + ylocal, elem.Klocal(xlocal, ylocal));
        if (nonzeroC != 0)
          C.write_value(3 * DIM * eIdx + xlocal, 3 * DIM * eIdx + ylocal, elem.Clocal(xlocal, ylocal));
      }
    }

    // (0a)
    // Initialize vector F
    F[3 * DIM * eIdx + 0] = elem.Flocal[0];
    F[3 * DIM * eIdx + 1] = elem.Flocal[1];
    F[3 * DIM * eIdx + 2] = elem.Flocal[2];
    F[3 * DIM * eIdx + 3] = elem.Flocal[3];
    F[3 * DIM * eIdx + 4] = elem.Flocal[4];
    F[3 * DIM * eIdx + 5] = elem.Flocal[5];
  }

  // ----------------------------------------------------------------------------------------------------------------

  int endnt = static_cast<int>(endtime / dt);
  int nt = 1;
  while (nt <= endnt) {
    float t = nt*dt;
    std::cout << "======= Time iteration #" << nt << "/" << endnt << " =======" << std::endl;
    std::cout << "=========== Time " << t << " ===========" << std::endl << std::endl;
    x_pred = x;
    x_pred.add_weighted(vel, 1.0f, dt);
    x_pred.add_weighted(res, 1.0f, 0.5f * (1.0f - beta2) * dt*dt);

    vel_pred = vel;
    vel_pred.add_weighted(res, 1.0f, (1.0f - beta1) * dt);

    b = C.MultiplyByVector(vel_pred, grid_size);

    b.add(K.MultiplyByVector(x_pred, grid_size));

    b.scale(-1.0f);

    // update F
    for (std::vector<BoundaryEdge>::iterator it = FEMdata.boundary.begin(); it != FEMdata.boundary.end(); ++it) {
      int eIdx = it->adj_elem1;
      if (DIM == 2)
        FEMdata.elements[eIdx].CalculateFlocal2D(*it, FEMdata.nodes, t);
      else if (DIM == 3)
        FEMdata.elements[eIdx].CalculateFlocal3D(*it, FEMdata.nodes, t);
            F[3 * DIM * eIdx + 0] = FEMdata.elements[eIdx].Flocal[0];
            F[3 * DIM * eIdx + 1] = FEMdata.elements[eIdx].Flocal[1];
            F[3 * DIM * eIdx + 2] = FEMdata.elements[eIdx].Flocal[2];
            F[3 * DIM * eIdx + 3] = FEMdata.elements[eIdx].Flocal[3];
            F[3 * DIM * eIdx + 4] = FEMdata.elements[eIdx].Flocal[4];
            F[3 * DIM * eIdx + 5] = FEMdata.elements[eIdx].Flocal[5];
    }

    b.add(F);  // Seems to be the same (that is, F), see formula (17.15)

        // Think how not to use loadVectors! Too dificult!
        std::unordered_map <int, MyArray> loadVectors;
        loadVectors.clear();    // in order to GetMapElement2Loadvector, because if not cleared the values are added
        // instead of assigned. See if-statement in for-loop in the function's body
        GetMapElement2Loadvector(FEMdata, loadVectors, t);
        std::cout << "loads\n";
        for (auto& it: loadVectors) {
          b[3 * DIM * it.first + 0] += loadVectors[it.first][0];
          b[3 * DIM * it.first + 1] += loadVectors[it.first][1];
          b[3 * DIM * it.first + 2] += loadVectors[it.first][2];
          b[3 * DIM * it.first + 3] += loadVectors[it.first][3];
          b[3 * DIM * it.first + 4] += loadVectors[it.first][4];
          b[3 * DIM * it.first + 5] += loadVectors[it.first][5];
          loadVectors[it.first].ShowNonzero();
          std::cout << "\n";
        }

    std::unordered_map <int, std::vector<int>> nodeAdjElem;
    CalculateNodeAdjElem(FEMdata, nodeAdjElem);

    ApplyConstraints_EbE(FEMdata);

    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
        it->Alocal = it->Mlocal.weightedSum(it->Clocal, 1.0f, beta1 * dt);
        it->Alocal = it->Mlocal;        // ????
        it->Alocal = it->Alocal.weightedSum(it->Klocal, 1.0f, 0.5f * beta2 * dt*dt);
    }

    AssignLoadElement(FEMdata, nodeAdjElem);

    //std::unordered_map <int, MyArray> loadVectors;

    // Set initial condition. Assumed that zero. ToDO: add initial conditions to prepared_meshes

    // ----------------------------------------------------------------------------------------------------------------
    int n_elems  = FEMdata.elementsCount;
    int n_gl_dofs = FEMdata.nodesCount * DIM;

    int grid_size = 3 * DIM * n_elems;
    MyArray x(grid_size), x_pred(grid_size), vel(grid_size),
            vel_pred(grid_size), b(grid_size), res(grid_size),
            F(grid_size);
    MyArray n_adjelem(n_gl_dofs);
    SparseMatrixCOO Q (grid_size);

    int nonzeroK = 0, nonzeroA = 0, nonzeroC = 0;
    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
        nonzeroK += it->Klocal.CountNonzero();
        nonzeroA += it->Alocal.CountNonzero();
        nonzeroC += it->Clocal.CountNonzero();
    }

    for (int eIdx = 0; eIdx < n_elems; ++eIdx) {
      FEMdata.elements[eIdx].blocal[0] = b[3 * DIM * eIdx + 0];
      FEMdata.elements[eIdx].blocal[1] = b[3 * DIM * eIdx + 1];
      FEMdata.elements[eIdx].blocal[2] = b[3 * DIM * eIdx + 2];
      FEMdata.elements[eIdx].blocal[3] = b[3 * DIM * eIdx + 3];
      FEMdata.elements[eIdx].blocal[4] = b[3 * DIM * eIdx + 4];
      FEMdata.elements[eIdx].blocal[5] = b[3 * DIM * eIdx + 5];
    }

    if (beta2 == 0.0f) {        // explicit scheme
      solve_diag_vec(FEMdata, res, false);
    } else {                      // implicit scheme
      PCG_EbE_vec_cpu(FEMdata, res, false, 1e-10f);
    }

    x = x_pred;
    x.add_weighted(res, 1.0f, 0.5f * beta2 * dt*dt);
    vel = vel_pred;
    vel.add_weighted(res, 1.0f, beta1 * dt);

    ++nt;
    std::cout << std::endl;
  }

  FEMdata.displacements = Q.MultiplyTransposedByVector(x, n_gl_dofs);
  FEMdata.displacements = FEMdata.displacements.divideByElementwise(n_adjelem);
}

void CalculateFEM_dyn_relaxation(FEMdataKeeper &FEMdata, float rho, float damping_alpha, float damping_beta, float dt, float beta1, float beta2, float eps) {
    // Zienkiewicz, Taylor, Zhu "The Finite Element Method: Its Basis and Fundamentals" 6th edition 17.3.3 GN22 (page 608)
    CheckRunTime(__func__)
    //const float beta1 = 0.55f; //0.55f;
    //const float beta2 = 0.5f; //0.5f;
    // beta2 = 0.0 -- explicit scheme (assuming both M and C are diagonal -- make sure to lump mass matrix!)
    // implicit scheme: beta1 >= beta2 >= 1/2
    int DIM = FEMdata.DIM;
    assert(beta2 == 0.0f || (beta1 >= beta2 && beta2 >= 0.5f));
    bool isLumping = (beta2 == 0.0f);

    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
        it->CalculateKlocal(FEMdata.D, FEMdata.nodes);
        it->CalculateMlocal(rho, FEMdata.nodes, isLumping);
        it->CalculateClocal(damping_alpha, damping_beta);
    }

    for (std::vector<BoundaryEdge>::iterator it = FEMdata.boundary.begin(); it != FEMdata.boundary.end(); ++it) {
      if (DIM == 2)
        FEMdata.elements[it->adj_elem1].CalculateFlocal2D(*it, FEMdata.nodes, 0.0f);
      else if (DIM == 3)
        FEMdata.elements[it->adj_elem1].CalculateFlocal3D(*it, FEMdata.nodes, 0.0f);
    }

    std::unordered_map <int, std::vector<int>> nodeAdjElem;
    CalculateNodeAdjElem(FEMdata, nodeAdjElem);

    ApplyConstraints_EbE(FEMdata);

    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
        it->Alocal = it->Mlocal.weightedSum(it->Clocal, 1.0f, beta1 * dt);
        it->Alocal = it->Mlocal;        // ????
        it->Alocal = it->Alocal.weightedSum(it->Klocal, 1.0f, 0.5f * beta2 * dt*dt);
    }

    AssignLoadElement(FEMdata, nodeAdjElem);

    //std::unordered_map <int, MyArray> loadVectors;

    // Set initial condition. Assumed that zero. ToDO: add initial conditions to prepared_meshes

    // ----------------------------------------------------------------------------------------------------------------
    int n_elems  = FEMdata.elementsCount;
    int n_gl_dofs = FEMdata.nodesCount * DIM;

    int grid_size = 3 * DIM * n_elems;
    MyArray x(grid_size), x_pred(grid_size), vel(grid_size),
            vel_pred(grid_size), b(grid_size), res(grid_size),
            F(grid_size);
    MyArray n_adjelem(n_gl_dofs);
    SparseMatrixCOO Q (grid_size);

    int nonzeroK = 0, nonzeroA = 0, nonzeroC = 0;
    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
        nonzeroK += it->Klocal.CountNonzero();
        nonzeroA += it->Alocal.CountNonzero();
        nonzeroC += it->Clocal.CountNonzero();
    }

    FEMdata.nonzeroMatrixNumCount = nonzeroA;

    SparseMatrixCOO K(nonzeroK), C(nonzeroC);
    for (int eIdx = 0; eIdx < n_elems; ++eIdx) {
        Element elem = FEMdata.elements[eIdx];

        Q.write_value(3 * DIM * eIdx + 0 * DIM + 0, elem.nodesIds[0] * DIM + 0, 1);
        Q.write_value(3 * DIM * eIdx + 0 * DIM + 1, elem.nodesIds[0] * DIM + 1, 1);
        Q.write_value(3 * DIM * eIdx + 1 * DIM + 0, elem.nodesIds[1] * DIM + 0, 1);
        Q.write_value(3 * DIM * eIdx + 1 * DIM + 1, elem.nodesIds[1] * DIM + 1, 1);
        Q.write_value(3 * DIM * eIdx + 2 * DIM + 0, elem.nodesIds[2] * DIM + 0, 1);
        Q.write_value(3 * DIM * eIdx + 2 * DIM + 1, elem.nodesIds[2] * DIM + 1, 1);

        n_adjelem[elem.nodesIds[0] * DIM + 0] += 1.0f;
        n_adjelem[elem.nodesIds[0] * DIM + 1] += 1.0f;
        n_adjelem[elem.nodesIds[1] * DIM + 0] += 1.0f;
        n_adjelem[elem.nodesIds[1] * DIM + 1] += 1.0f;
        n_adjelem[elem.nodesIds[2] * DIM + 0] += 1.0f;
        n_adjelem[elem.nodesIds[2] * DIM + 1] += 1.0f;

        for (int xlocal = 0; xlocal < 3 * DIM; ++xlocal) {
            for (int ylocal = 0; ylocal < 3 * DIM; ++ylocal) {
                K.write_value(3 * DIM * eIdx + xlocal, 3 * DIM * eIdx + ylocal, elem.Klocal(xlocal, ylocal));
                if (nonzeroC != 0)
                    C.write_value(3 * DIM * eIdx + xlocal, 3 * DIM * eIdx + ylocal, elem.Clocal(xlocal, ylocal));
            }
        }

        // (0a)
        // Initialize vector F
        F[3 * DIM * eIdx + 0] = elem.Flocal[0];
        F[3 * DIM * eIdx + 1] = elem.Flocal[1];
        F[3 * DIM * eIdx + 2] = elem.Flocal[2];
        F[3 * DIM * eIdx + 3] = elem.Flocal[3];
        F[3 * DIM * eIdx + 4] = elem.Flocal[4];
        F[3 * DIM * eIdx + 5] = elem.Flocal[5];
    }

    // ----------------------------------------------------------------------------------------------------------------

    int nt = 1;
    float cnorm_res, cnorm_vel;
    do {
        float t = nt*dt;
        std::cout << "======= Time iteration #" << nt << " =======" << std::endl;
        std::cout << "=========== Time " << t << " ===========" << std::endl << std::endl;
        x_pred = x;
        x_pred.add_weighted(vel, 1.0f, dt);
        x_pred.add_weighted(res, 1.0f, 0.5f * (1.0f - beta2) * dt*dt);

        vel_pred = vel;
        vel_pred.add_weighted(res, 1.0f, (1.0f - beta1) * dt);

        b = C.MultiplyByVector(vel_pred, grid_size);

        b.add(K.MultiplyByVector(x_pred, grid_size));

        b.scale(-1.0f);

        // update F
        for (std::vector<BoundaryEdge>::iterator it = FEMdata.boundary.begin(); it != FEMdata.boundary.end(); ++it) {
            int eIdx = it->adj_elem1;
            if (DIM == 2)
              FEMdata.elements[it->adj_elem1].CalculateFlocal2D(*it, FEMdata.nodes, t);
            else if (DIM == 3)
              FEMdata.elements[it->adj_elem1].CalculateFlocal3D(*it, FEMdata.nodes, t);
            F[3 * DIM * eIdx + 0] = FEMdata.elements[eIdx].Flocal[0];
            F[3 * DIM * eIdx + 1] = FEMdata.elements[eIdx].Flocal[1];
            F[3 * DIM * eIdx + 2] = FEMdata.elements[eIdx].Flocal[2];
            F[3 * DIM * eIdx + 3] = FEMdata.elements[eIdx].Flocal[3];
            F[3 * DIM * eIdx + 4] = FEMdata.elements[eIdx].Flocal[4];
            F[3 * DIM * eIdx + 5] = FEMdata.elements[eIdx].Flocal[5];
        }
        b.add(F);  // Seems to be the same (that is, F), see formula (17.15)

        // Think how not to use loadVectors! Too dificult!
        std::unordered_map <int, MyArray> loadVectors;
        loadVectors.clear();    // in order to GetMapElement2Loadvector, because if not cleared the values are added
                                // instead of assigned. See if-statement in for-loop in the function's body
        GetMapElement2Loadvector(FEMdata, loadVectors, t);
        std::cout << "loads\n";
        for (auto& it: loadVectors) {
            b[3 * DIM * it.first + 0] += loadVectors[it.first][0];
            b[3 * DIM * it.first + 1] += loadVectors[it.first][1];
            b[3 * DIM * it.first + 2] += loadVectors[it.first][2];
            b[3 * DIM * it.first + 3] += loadVectors[it.first][3];
            b[3 * DIM * it.first + 4] += loadVectors[it.first][4];
            b[3 * DIM * it.first + 5] += loadVectors[it.first][5];
            loadVectors[it.first].ShowNonzero();
            std::cout << "\n";
        }

        for (int eIdx = 0; eIdx < n_elems; ++eIdx) {
            FEMdata.elements[eIdx].blocal[0] = b[3 * DIM * eIdx + 0];
            FEMdata.elements[eIdx].blocal[1] = b[3 * DIM * eIdx + 1];
            FEMdata.elements[eIdx].blocal[2] = b[3 * DIM * eIdx + 2];
            FEMdata.elements[eIdx].blocal[3] = b[3 * DIM * eIdx + 3];
            FEMdata.elements[eIdx].blocal[4] = b[3 * DIM * eIdx + 4];
            FEMdata.elements[eIdx].blocal[5] = b[3 * DIM * eIdx + 5];
        }

        if (beta2 == 0.0f) {        // explicit scheme
            solve_diag_vec(FEMdata, res, false);
        } else {                      // implicit scheme
            PCG_EbE_vec_cpu(FEMdata, res, false, 1e-10f);
        }

        x = x_pred;
        x.add_weighted(res, 1.0f, 0.5f * beta2 * dt*dt);
        vel = vel_pred;
        vel.add_weighted(res, 1.0f, beta1 * dt);

        ++nt;
        cnorm_res = res.CNorm();
        cnorm_vel = vel.CNorm();
        std::cout << "C-norm res = " << cnorm_res << "\nC-norm vel = " << cnorm_vel << std::endl;
        std::cout << std::endl;
        if (nt >= 10000) break;
    } while (!((cnorm_res < eps) && (cnorm_vel < eps)));

    FEMdata.displacements = Q.MultiplyTransposedByVector(x, n_gl_dofs);
    FEMdata.displacements = FEMdata.displacements.divideByElementwise(n_adjelem);
}

void CalculateFEM_dyn_explicit(FEMdataKeeper &FEMdata, float rho, float damping_alpha, float damping_beta, float endtime, float dt, float beta1) {
    // Zienkiewicz, Taylor, Zhu "The Finite Element Method: Its Basis and Fundamentals" 6th edition 17.3.3 GN22 (page 608)
    CheckRunTime(__func__)

    int DIM = FEMdata.DIM;

    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
        it->CalculateKlocal(FEMdata.D, FEMdata.nodes);
        it->CalculateMlocal(rho, FEMdata.nodes, true); // Lumping on
        it->CalculateClocal(damping_alpha, damping_beta);
    }

    for (std::vector<BoundaryEdge>::iterator it = FEMdata.boundary.begin(); it != FEMdata.boundary.end(); ++it) {
      if (DIM == 2)
        FEMdata.elements[it->adj_elem1].CalculateFlocal2D(*it, FEMdata.nodes, 0.0f);
      else if (DIM == 3)
        FEMdata.elements[it->adj_elem1].CalculateFlocal3D(*it, FEMdata.nodes, 0.0f);
    }

    std::unordered_map <int, std::vector<int>> nodeAdjElem;
    CalculateNodeAdjElem(FEMdata, nodeAdjElem);

    ApplyConstraints_EbE(FEMdata);

    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
        it->Alocal = it->Mlocal.weightedSum(it->Clocal, 1.0f, beta1 * dt);
        it->Alocal = it->Mlocal;        // ???
    }

    AssignLoadElement(FEMdata, nodeAdjElem);

    //std::unordered_map <int, MyArray> loadVectors;

    // Set initial condition. Assumed that zero. ToDO: add initial conditions to prepared_meshes

    // ----------------------------------------------------------------------------------------------------------------
    int n_elems  = FEMdata.elementsCount;
    int n_gl_dofs = FEMdata.nodesCount * DIM;

    int grid_size = 3 * DIM * n_elems;
    MyArray x(grid_size), vel(grid_size),
            b(grid_size), res(grid_size),
            F(grid_size);
    MyArray n_adjelem(n_gl_dofs);
    SparseMatrixCOO Q (grid_size);

    int nonzeroK = 0, nonzeroA = 0, nonzeroC = 0;
    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
        nonzeroK += it->Klocal.CountNonzero();
        nonzeroA += it->Alocal.CountNonzero();
        nonzeroC += it->Clocal.CountNonzero();
    }

    FEMdata.nonzeroMatrixNumCount = nonzeroA;

    SparseMatrixCOO K(nonzeroK), C(nonzeroC);
    for (int eIdx = 0; eIdx < n_elems; ++eIdx) {
        Element elem = FEMdata.elements[eIdx];

        Q.write_value(3 * DIM * eIdx + 0 * DIM + 0, elem.nodesIds[0] * DIM + 0, 1);
        Q.write_value(3 * DIM * eIdx + 0 * DIM + 1, elem.nodesIds[0] * DIM + 1, 1);
        Q.write_value(3 * DIM * eIdx + 1 * DIM + 0, elem.nodesIds[1] * DIM + 0, 1);
        Q.write_value(3 * DIM * eIdx + 1 * DIM + 1, elem.nodesIds[1] * DIM + 1, 1);
        Q.write_value(3 * DIM * eIdx + 2 * DIM + 0, elem.nodesIds[2] * DIM + 0, 1);
        Q.write_value(3 * DIM * eIdx + 2 * DIM + 1, elem.nodesIds[2] * DIM + 1, 1);

        n_adjelem[elem.nodesIds[0] * DIM + 0] += 1.0f;
        n_adjelem[elem.nodesIds[0] * DIM + 1] += 1.0f;
        n_adjelem[elem.nodesIds[1] * DIM + 0] += 1.0f;
        n_adjelem[elem.nodesIds[1] * DIM + 1] += 1.0f;
        n_adjelem[elem.nodesIds[2] * DIM + 0] += 1.0f;
        n_adjelem[elem.nodesIds[2] * DIM + 1] += 1.0f;

        for (int xlocal = 0; xlocal < 3 * DIM; ++xlocal) {
            for (int ylocal = 0; ylocal < 3 * DIM; ++ylocal) {
                K.write_value(3 * DIM * eIdx + xlocal, 3 * DIM * eIdx + ylocal, elem.Klocal(xlocal, ylocal));
                if (nonzeroC != 0)
                    C.write_value(3 * DIM * eIdx + xlocal, 3 * DIM * eIdx + ylocal, elem.Clocal(xlocal, ylocal));
            }
        }

        // (0a)
        // Initialize vector F
        F[3 * DIM * eIdx + 0] = elem.Flocal[0];
        F[3 * DIM * eIdx + 1] = elem.Flocal[1];
        F[3 * DIM * eIdx + 2] = elem.Flocal[2];
        F[3 * DIM * eIdx + 3] = elem.Flocal[3];
        F[3 * DIM * eIdx + 4] = elem.Flocal[4];
        F[3 * DIM * eIdx + 5] = elem.Flocal[5];
    }

    // ----------------------------------------------------------------------------------------------------------------

    int endnt = static_cast<int>(endtime / dt);
    int nt = 1;
    while (nt <= endnt) {
        float t = nt*dt;
        std::cout << "======= Time iteration #" << nt << "/" << endnt << " =======" << std::endl;
        std::cout << "=========== Time " << t << " ===========" << std::endl << std::endl;
        x.add_weighted(vel, 1.0f, dt);
        x.add_weighted(res, 1.0f, 0.5f * dt*dt);

        vel.add_weighted(res, 1.0f, (1.0f - beta1) * dt);

        b = C.MultiplyByVector(vel, grid_size);

        b.add(K.MultiplyByVector(x, grid_size));

        b.scale(-1.0f);

        // update F
        for (std::vector<BoundaryEdge>::iterator it = FEMdata.boundary.begin(); it != FEMdata.boundary.end(); ++it) {
            int eIdx = it->adj_elem1;
            if (DIM == 2)
              FEMdata.elements[it->adj_elem1].CalculateFlocal2D(*it, FEMdata.nodes, t);
            else if (DIM == 3)
              FEMdata.elements[it->adj_elem1].CalculateFlocal3D(*it, FEMdata.nodes, t);
            F[3 * DIM * eIdx + 0] = FEMdata.elements[eIdx].Flocal[0];
            F[3 * DIM * eIdx + 1] = FEMdata.elements[eIdx].Flocal[1];
            F[3 * DIM * eIdx + 2] = FEMdata.elements[eIdx].Flocal[2];
            F[3 * DIM * eIdx + 3] = FEMdata.elements[eIdx].Flocal[3];
            F[3 * DIM * eIdx + 4] = FEMdata.elements[eIdx].Flocal[4];
            F[3 * DIM * eIdx + 5] = FEMdata.elements[eIdx].Flocal[5];
        }
        b.add(F);  // Seems to be the same (that is, F), see formula (17.15)

        // Think how not to use loadVectors! Too dificult!
        std::unordered_map <int, MyArray> loadVectors;
        loadVectors.clear();    // in order to GetMapElement2Loadvector, because if not cleared the values are added
                                // instead of assigned. See if-statement in for-loop in the function's body
        GetMapElement2Loadvector(FEMdata, loadVectors, t);
        std::cout << "loads\n";
        for (auto& it: loadVectors) {
            b[3 * DIM * it.first + 0] += loadVectors[it.first][0];
            b[3 * DIM * it.first + 1] += loadVectors[it.first][1];
            b[3 * DIM * it.first + 2] += loadVectors[it.first][2];
            b[3 * DIM * it.first + 3] += loadVectors[it.first][3];
            b[3 * DIM * it.first + 4] += loadVectors[it.first][4];
            b[3 * DIM * it.first + 5] += loadVectors[it.first][5];
            loadVectors[it.first].ShowNonzero();
            std::cout << "\n";
        }

        for (int eIdx = 0; eIdx < n_elems; ++eIdx) {
            FEMdata.elements[eIdx].blocal[0] = b[3 * DIM * eIdx + 0];
            FEMdata.elements[eIdx].blocal[1] = b[3 * DIM * eIdx + 1];
            FEMdata.elements[eIdx].blocal[2] = b[3 * DIM * eIdx + 2];
            FEMdata.elements[eIdx].blocal[3] = b[3 * DIM * eIdx + 3];
            FEMdata.elements[eIdx].blocal[4] = b[3 * DIM * eIdx + 4];
            FEMdata.elements[eIdx].blocal[5] = b[3 * DIM * eIdx + 5];
        }

        solve_diag_vec(FEMdata, res, false);

        vel.add_weighted(res, 1.0f, beta1 * dt);

        ++nt;
        std::cout << std::endl;
    }

    FEMdata.displacements = Q.MultiplyTransposedByVector(x, n_gl_dofs);
    FEMdata.displacements = FEMdata.displacements.divideByElementwise(n_adjelem);
}

//float gpuDotProduct(float *a, float *b, int size) {
//    thrust::device_vector<float> d_a(a, a + size);
//    thrust::device_vector<float> d_b(b, b + size);
//    return thrust::inner_product(thrust::cuda::par, d_a.begin(), d_a.end(), d_b.begin(), 0.0f);
//}

// 2D only
void SmoothResults(std::string stress_component, MyArray &SmoothStress, std::vector<MyArray> Stress,
                   int nodesCount, std::vector<MyArray> &nodes, std::vector<Element> elements) {
  Matrix C(nodesCount);
  MyArray R(nodesCount);

  int StressIdx = 0;
  if (stress_component == "xx") {
    StressIdx = 0;
  } else if (stress_component == "yy") {
    StressIdx = 1;
  } else if (stress_component == "xy") {
    StressIdx = 2;
  }

  float dlt, x1, y1, x2, y2, x3, y3;
  for (int i = 0; i < elements.size(); ++i) {
    x1 = nodes[0][elements[i].nodesIds[0]]; y1 = nodes[1][elements[i].nodesIds[0]];
    x2 = nodes[0][elements[i].nodesIds[1]]; y2 = nodes[1][elements[i].nodesIds[1]];
    x3 = nodes[0][elements[i].nodesIds[2]]; y3 = nodes[1][elements[i].nodesIds[2]];

    float a1 = x2 * y3 - x3 * y2;
    float b1 = y2 - y3;
    float c1 = x3 - x2;

    float a2 = x3 * y1 - x1 * y3;
    float b2 = y3 - y1;
    float c2 = x1 - x3;

    float a3 = x1 * y2 - x2 * y1;
    float b3 = y1 - y2;
    float c3 = x2 - x1;

    float x_c = (x1 + x2 + x3) / 3.0f;
    float y_c = (y1 + y2 + y3) / 3.0f;

    dlt = std::abs((x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3)) / 2.0f;

    //        R[elements[i].nodesIds[0]] += (a1 + b1 * x_c + c1 * y_c) * 0.5 * Stress[i][StressIdx];
    //        R[elements[i].nodesIds[1]] += (a2 + b2 * x_c + c2 * y_c) * 0.5 * Stress[i][StressIdx];
    //        R[elements[i].nodesIds[2]] += (a3 + b3 * x_c + c3 * y_c) * 0.5 * Stress[i][StressIdx];

    R[elements[i].nodesIds[0]] += Stress[i][StressIdx] * dlt * 3.0f;
    R[elements[i].nodesIds[1]] += Stress[i][StressIdx] * dlt * 3.0f;
    R[elements[i].nodesIds[2]] += Stress[i][StressIdx] * dlt * 3.0f;

    //        C(elements[i].nodesIds[0], elements[i].nodesIds[0]) += (a1*a1+(a1*b1+a1*b1)*x_c+(a1*c1+a1*c1)*y_c+b1*b1*x_c*x_c+c1*c1*y_c*y_c+(b1*c1+b1*c1)*x_c*y_c) / (4 * dlt);
    //        C(elements[i].nodesIds[0], elements[i].nodesIds[1]) += (a1*a2+(a1*b2+a2*b1)*x_c+(a1*c2+a2*c1)*y_c+b1*b2*x_c*x_c+c1*c2*y_c*y_c+(b1*c2+b2*c1)*x_c*y_c) / (4 * dlt);
    //        C(elements[i].nodesIds[0], elements[i].nodesIds[2]) += (a1*a3+(a1*b3+a3*b1)*x_c+(a1*c3+a3*c1)*y_c+b1*b3*x_c*x_c+c1*c3*y_c*y_c+(b1*c3+b3*c1)*x_c*y_c) / (4 * dlt);
    //        C(elements[i].nodesIds[1], elements[i].nodesIds[0]) += (a2*a1+(a2*b1+a1*b2)*x_c+(a2*c1+a1*c2)*y_c+b2*b1*x_c*x_c+c2*c1*y_c*y_c+(b2*c1+b1*c2)*x_c*y_c) / (4 * dlt);
    //        C(elements[i].nodesIds[1], elements[i].nodesIds[1]) += (a2*a2+(a2*b2+a2*b2)*x_c+(a2*c2+a2*c2)*y_c+b2*b2*x_c*x_c+c2*c2*y_c*y_c+(b2*c2+b2*c2)*x_c*y_c) / (4 * dlt);
    //        C(elements[i].nodesIds[1], elements[i].nodesIds[2]) += (a2*a3+(a2*b3+a3*b2)*x_c+(a2*c3+a3*c2)*y_c+b2*b3*x_c*x_c+c2*c3*y_c*y_c+(b2*c3+b3*c2)*x_c*y_c) / (4 * dlt);
    //        C(elements[i].nodesIds[2], elements[i].nodesIds[0]) += (a3*a1+(a3*b1+a1*b3)*x_c+(a3*c1+a1*c3)*y_c+b3*b1*x_c*x_c+c3*c1*y_c*y_c+(b3*c1+b1*c3)*x_c*y_c) / (4 * dlt);
    //        C(elements[i].nodesIds[2], elements[i].nodesIds[1]) += (a3*a2+(a3*b2+a2*b3)*x_c+(a3*c2+a2*c3)*y_c+b3*b2*x_c*x_c+c3*c2*y_c*y_c+(b3*c2+b2*c3)*x_c*y_c) / (4 * dlt);
    //        C(elements[i].nodesIds[2], elements[i].nodesIds[2]) += (a3*a3+(a3*b3+a3*b3)*x_c+(a3*c3+a3*c3)*y_c+b3*b3*x_c*x_c+c3*c3*y_c*y_c+(b3*c3+b3*c3)*x_c*y_c) / (4 * dlt);

    C(elements[i].nodesIds[0], elements[i].nodesIds[0]) += dlt;
    C(elements[i].nodesIds[0], elements[i].nodesIds[1]) += dlt;
    C(elements[i].nodesIds[0], elements[i].nodesIds[2]) += dlt;
    C(elements[i].nodesIds[1], elements[i].nodesIds[0]) += dlt;
    C(elements[i].nodesIds[1], elements[i].nodesIds[1]) += dlt;
    C(elements[i].nodesIds[1], elements[i].nodesIds[2]) += dlt;
    C(elements[i].nodesIds[2], elements[i].nodesIds[0]) += dlt;
    C(elements[i].nodesIds[2], elements[i].nodesIds[1]) += dlt;
    C(elements[i].nodesIds[2], elements[i].nodesIds[2]) += dlt;
  }
  C.LU_solve(R, SmoothStress, nodesCount);
}

void MakeResults(FEMdataKeeper FEMdata, ResultsDataKeeper &RESdata) {
  //POSTPROCESSING

  CalculateStressAndDeformation(FEMdata.DIM, RESdata.Deformation,
                                RESdata.Stress,
                                RESdata.epsilon_mises,
                                RESdata.sigma_mises,
                                FEMdata.D, FEMdata.elements, FEMdata.displacements, FEMdata.all_B);

  float fixed_value = 1.0;// x -3.0; y -3.0
  float a = 0.0f;// x -3.0 y -4.0
  float b = 8.0f;// x 5.0 y -3.0;

  if (RESdata.withStressAlongAxis) {
    CalculateStressAlongAxis(RESdata.StressComponents,
                           "x", "xy", fixed_value, a, b,
                           RESdata.Stress, FEMdata.nodes, FEMdata.elements);
  }
  if (RESdata.withSmooth) {
    SmoothResults("xx", RESdata.SmoothStress, RESdata.Stress, FEMdata.nodesCount, FEMdata.nodes, FEMdata.elements);
    CalculateStressAlongAxisSmooth(RESdata.StressComponentsSmooth,
                                   "x", fixed_value, a, b, RESdata.SmoothStress, FEMdata.nodes, FEMdata.elements);
  }

  a = 0.0f;
  b = 6.0f;
  float k = -0.6655f;
  float m = -0.0035f;

  if (RESdata.withMises) {
    CalculateMisesAlongLineMises(RESdata.MisesComponents,
                                 k, m, a, b, RESdata.sigma_mises, FEMdata.nodes, FEMdata.elements);
  }
}

void WriteResults(FEMdataKeeper FEMdata, ResultsDataKeeper RESdata, std::string output_vtk, bool PRINT_DEBUG_INFO) {
  std::string path_stress = FEMdata.res_dir + "/output/out_stress_" + FEMdata.get_name() + ".txt";
  if (RESdata.withStressAlongAxis) {
    std::cout << "StressComponents Size = " << RESdata.StressComponents.size() << "\n";
    fstream out1;

    out1.open(path_stress, fstream::out);
    out1 << FEMdata.nodesCount << " " << FEMdata.elementsCount << "\n";
    for (int i = 0; i < RESdata.StressComponents.size(); i+=2) {
      out1 << RESdata.StressComponents[i] << " " << RESdata.StressComponents[i + 1] << "\n";
    }
  }

  if (RESdata.withSmooth) {
    fstream out2;
    std::string path_stress_smooth = FEMdata.res_dir + "/output/out_stress_" + FEMdata.get_name() + "_smooth.txt";
    if (PRINT_DEBUG_INFO) {
      std::cout << "StressComponentsSmooth Size = " << RESdata.StressComponentsSmooth.size() << "\n";
    }
    out2.open(path_stress_smooth, fstream::out);
    out2 << FEMdata.nodesCount << " " << FEMdata.elementsCount << "\n";
    for (int i = 0; i < RESdata.StressComponentsSmooth.size(); i+=2) {
      out2 << RESdata.StressComponentsSmooth[i] << " " << RESdata.StressComponentsSmooth[i + 1] << "\n";
    }
  }

  if (RESdata.withMises) {
    fstream out3;
    std::string path_stress_mises = FEMdata.res_dir + "/output/out_stress_" + FEMdata.get_name() + "_mises.txt";
    if (PRINT_DEBUG_INFO) {
      std::cout << "MisesComponents Size = " << RESdata.MisesComponents.size() << "\n";
    }
    out3.open(path_stress_mises, fstream::out);
    out3 << FEMdata.nodesCount << " " << FEMdata.elementsCount << "\n";
    for (int i = 0; i < RESdata.MisesComponents.size(); i+=2) {
      out3 << RESdata.MisesComponents[i] << " " << RESdata.MisesComponents[i + 1] << "\n";
    }
  }

  if (FEMdata.DIM == 2) {
    MakeVTKfile2D(output_vtk, FEMdata.nodes, FEMdata.elements,
                FEMdata.displacements, RESdata.Stress, RESdata.sigma_mises, RESdata.Deformation, RESdata.epsilon_mises, RESdata.SmoothStress);
  } else if (FEMdata.DIM == 3) {
    MakeVTKfile3D(output_vtk, FEMdata.nodes, FEMdata.elements,
                FEMdata.displacements, RESdata.Stress, RESdata.sigma_mises, RESdata.Deformation, RESdata.epsilon_mises);
  }
}
