#include "femfunc.h"

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
  const int DIM = FEMdata.DIM;
  for (std::vector<Load>::iterator it = FEMdata.loads.begin(); it != FEMdata.loads.end(); ++it) {
    assert(it->elem != -1);
    const int loadGlobalNode  = it->dof / DIM;
    const int loadFreedomAxis = it->dof % DIM;
    Element *elem = &FEMdata.elements[it->elem];
    elem->Flocal[elem->Global2LocalNode(loadGlobalNode) * DIM + loadFreedomAxis] += it->value;
  }
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

void CalculateNodeAdjElem(FEMdataKeeper &FEMdata, std::unordered_map <int, std::vector<int>> &a) {
  for (int n = 0; n < FEMdata.elements.size(); ++n) {
    a[FEMdata.elements[n].nodesIds[0]].push_back(n);
    a[FEMdata.elements[n].nodesIds[1]].push_back(n);
    a[FEMdata.elements[n].nodesIds[2]].push_back(n);
    if (FEMdata.DIM == 3) {
      a[FEMdata.elements[n].nodesIds[3]].push_back(n);
    }
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

void gpuCalculateFEM_EbE_vec(FEMdataKeeper &FEMdata, bool PRINT_DEBUG_INFO) {
  CheckRunTime(__func__)

  // Note that local matrices are calculated on GPU, see function gpuPCG_EbE_vec

  // ToDO: Unite CalculateFlocal2D and CalculateFlocal3D
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

  ApplyConstraints_EbE(FEMdata);

  gpuPCG_EbE_vec(FEMdata, FEMdata.displacements, true, 1e-4f, PRINT_DEBUG_INFO);
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
  // Zienkiewicz, Taylor, Zhu "The Finite Element Method: Its Basis and Fundamentals" 6th edition 17.3.3 GN22 (page 608)
  CheckRunTime(__func__)
  //assert(beta2 >= beta1 && beta1 >= 0.5f);
  bool isExplicit = beta2 == 0;
  if (isExplicit)
    assert(damping_beta == 0.0f);

  int n_elems = FEMdata.elementsCount;
  int n_gl_dofs = FEMdata.nodesCount * FEMdata.DIM;
  int grid_size = 3 * FEMdata.DIM * n_elems;
  bool doAssemblyRes = false;
  bool isLumped = isExplicit;
  bool isDamping = !(damping_alpha == 0.0f && damping_beta == 0.0f);
  float eps_PCG   = 1e-4f;
  float eps_relax;
  bool is_relax;

  for (std::vector<BoundaryEdge>::iterator it = FEMdata.boundary.begin(); it != FEMdata.boundary.end(); ++it) {
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

  gpuDataKeeper_DYN_DAMP gpu_data(FEMdata.DIM, FEMdata.elementsCount, FEMdata.nodesCount, doAssemblyRes, isLumped, damping_alpha, damping_beta);
  if (isDamping) {
    if (isExplicit)
      gpu_data.set_SLAU_matrix_coefs(1.0f, 0.0f, beta1 * dt);
    else
      gpu_data.set_SLAU_matrix_coefs(1.0f, 0.5f * beta2 * dt*dt, beta1 * dt);
  } else if (!isExplicit) {
    gpu_data.set_SLAU_matrix_coefs(1.0f, 0.5f * beta2 * dt*dt, 0.0f);
  }
  copyElementsAndFlocals(FEMdata, gpu_data);
  gpuGenerateMask(gpu_data, FEMdata.DIM, FEMdata.elementsCount);
  gpuCountNAdjElem(gpu_data, grid_size);
  gpuCalculateKlocal2(gpu_data, FEMdata);
  gpu_data.copyBmatrixToHost(FEMdata.all_B.get_data());
  gpuCalculateMlocal(gpu_data, FEMdata, rho);
  if (isDamping)
    gpuCalculateDiag_DAMP(gpu_data, n_elems);
  else if (!isExplicit)
    gpuCalculateDiag(gpu_data, n_elems);

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
        if (!is_relax) cout << "/" << endnt << " =======";
        std::cout << " =======";
        std::cout << "\n=========== Time " << t << " ===========\n\n";
      }

      gpuAddWeighted2(gpu_data.get_displ(), gpu_data.get_vel(), 1.0f, dt, grid_size);
      gpuAddWeighted2(gpu_data.get_displ(), gpu_data.get_x(), 1.0f, 0.5f * (1.0f - beta2) *  dt*dt, grid_size);
      gpuAddWeighted2(gpu_data.get_vel(), gpu_data.get_x(), 1.0f, (1.0f - beta1) * dt, grid_size);
      gpu_data.setZeroVec();

      gpuMultiplyMatrixByVec(gpu_data.get_Klocals(), gpu_data.get_displ(), gpu_data.get_r(), FEMdata.elementsCount);
      if (isDamping) {
        gpuMultiplyClocalByVec(gpu_data, gpu_data.get_vel(), gpu_data.get_tmp(), FEMdata.elementsCount);
        gpuAdd(gpu_data.get_r(), gpu_data.get_tmp(), grid_size);
      }

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

      if (isDamping) {
        if (isExplicit)
          gpuSolveDiag(gpu_data.get_diag(),  gpu_data.get_r(),
                       gpu_data.get_x(), gpu_data.get_mask(),
                       gpu_data.get_n_adjelem(), grid_size, n_gl_dofs, doAssemblyRes);
        else
          gpuPCG_EbE_vec_DYN_DAMP(FEMdata, gpu_data, doAssemblyRes, eps_PCG, PRINT_DEBUG_INFO);
       } else {
        if (isExplicit)
          gpuSolveDiag(gpu_data.get_diagM(), gpu_data.get_r(),
                       gpu_data.get_x(), gpu_data.get_mask(),
                       gpu_data.get_n_adjelem(), grid_size, n_gl_dofs, doAssemblyRes);
        else
          gpuPCG_EbE_vec_DYN(FEMdata, gpu_data, doAssemblyRes, eps_PCG, PRINT_DEBUG_INFO);
      }

      gpuAddWeighted2(gpu_data.get_vel(), gpu_data.get_x(), 1.0f, beta1 * dt, grid_size);
      if (!isExplicit) gpuAddWeighted2(gpu_data.get_displ(), gpu_data.get_x(), 1.0f, 0.5f * beta2 * dt*dt, grid_size);


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

void MakeResults(FEMdataKeeper &FEMdata, ResultsDataKeeper &RESdata) {
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

void WriteResults(FEMdataKeeper &FEMdata, ResultsDataKeeper &RESdata, std::string output_vtk, bool PRINT_DEBUG_INFO) {
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
