
#include "femfunc.h"

using namespace std;


float SetConstraints(int i, int j, float v, int index) {
    if (i == index || j == index) {
        return i == j ? 1.0 : 0.0;
    } else {
        return v;
    }
}

void ApplyConstraints(SparseMatrixCOO& K, const std::vector<Constraint>& constraints, int n) {
    CheckRunTime(__func__)
    std::vector<int> indicesToConstraint;

    for (std::vector<Constraint>::const_iterator it = constraints.begin(); it != constraints.end(); ++it) {
        if (it->type & Constraint::UX) {
            indicesToConstraint.push_back(3 * it->node + 0);
        }
        if (it->type & Constraint::UY) {
            indicesToConstraint.push_back(3 * it->node + 1);
        }
        if (it->type & Constraint::UZ) {
            indicesToConstraint.push_back(3 * it->node + 2);
        }
    }


//    for (int i = 0; i < indicesToConstraint.size(); i++) {
//        cout << indicesToConstraint[i] << " ";
//    }

    int k = 0;
    int i = 0;
    int threshhold = K.get_size();
    while (i < threshhold) {
    //for (int i = 0; i < K.get_size(); i++) {
        for (int j = 0; j < indicesToConstraint.size(); j++) {
            if (K.get_x(i) == indicesToConstraint[j] || K.get_y(i) == indicesToConstraint[j]) {
                if (K.get_x(i) == K.get_y(i)) {
                    K.set_value(K.get_x(i), K.get_y(i), 1.0);
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
    //cout << "!!!!" << k << " " << indicesToConstraint.size() << " ";
}

void CalculateStressAndDeformation(std::vector<MyArray> &Deformation,
                                   std::vector<MyArray> &Stress,
                                   std::vector<float> &epsilon_mises,
                                   std::vector<float> &sigma_mises,
                                   Matrix &D,
                                   std::vector<Element> &elements,
                                   MyArray &displacements,
                                   MyArray &all_B) {
    CheckRunTime(__func__)
    MyArray StressVector(6);
    MyArray DeformationVector(6);
    MyArray delta(12);

    Matrix B(6, 12);

    for (int i = 0; i < elements.size(); ++i) {
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

        for (int k = 0; k < 6; ++k)
              for (int j = 0; j < 12; ++j)
                B(k, j) = all_B[j + 12 * k + 6 * 12 * i];
        DeformationVector = B.Product(delta);
        StressVector = D.Product(DeformationVector);

        float sigma = sqrt(StressVector[0] * StressVector[0] - StressVector[0]
                    * StressVector[1] + StressVector[1] * StressVector[1] + 3.0 * StressVector[2] * StressVector[2]);
        sigma_mises.push_back(sigma);

        float epsilon = sqrt(DeformationVector[0] * DeformationVector[0] - DeformationVector[0]
                * DeformationVector[1] + DeformationVector[1] * DeformationVector[1] + 3.0 * DeformationVector[2] * DeformationVector[2]);
        epsilon_mises.push_back(epsilon);

        Deformation.push_back(DeformationVector);
        Stress.push_back(StressVector);

        //DeformationVector.Show();
        //StressVector.Show();
        //cout << endl;
    }
}

//void CalculateFiniteElementMethod(FEMdataKeeper &FEMdata) {
//    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
//        it->CalculateKlocal(FEMdata.D, FEMdata.nodesX, FEMdata.nodesY, FEMdata.nodesZ);
//    }

//    SparseMatrixCOO globalK = AssemblyStiffnessMatrix(FEMdata);

//    ApplyConstraints(globalK, FEMdata.constraints, FEMdata.loads.get_size());

//    globalK.SortIt();

//    cout << "nonzero = " << globalK.CountNonZero() << endl;

//    globalK.CGM_solve(FEMdata.loads, FEMdata.displacements, FEMdata.loads.get_size(), 1e-10);
//}

void MakeResults(FEMdataKeeper &FEMdata, ResultsDataKeeper &RESdata) {
  //POSTPROCESSING

  CalculateStressAndDeformation(RESdata.Deformation,
                                RESdata.Stress,
                                RESdata.epsilon_mises,
                                RESdata.sigma_mises,
                                FEMdata.D, FEMdata.elements, FEMdata.displacements, FEMdata.all_B);

//  float fixed_value = 1.0;// x -3.0; y -3.0
//  float a = 0.0f;// x -3.0 y -4.0
//  float b = 8.0f;// x 5.0 y -3.0;

//  CalculateStressAlongAxis(RESdata.StressComponents,
//                           "x", "xy", fixed_value, a, b,
//                           RESdata.Stress, FEMdata.nodesX, FEMdata.nodesY, FEMdata.elements);

//  if (RESdata.withSmooth) {
//    SmoothResults("xx", RESdata.SmoothStress, RESdata.Stress, FEMdata.nodesCount, FEMdata.nodesX, FEMdata.nodesY, FEMdata.elements);
//    CalculateStressAlongAxisSmooth(RESdata.StressComponentsSmooth,
//                                   "x", fixed_value, a, b, RESdata.SmoothStress, FEMdata.nodesX, FEMdata.nodesY, FEMdata.elements);
//  }

//  a = 0.0f;
//  b = 6.0f;
//  float k = -0.6655f;
//  float m = -0.0035f;

//  if (RESdata.withMises) {
//    CalculateMisesAlongLineMises(RESdata.MisesComponents,
//                                 k, m, a, b, RESdata.sigma_mises, FEMdata.nodesX, FEMdata.nodesY, FEMdata.elements);
//  }
}

void WriteResults(FEMdataKeeper &FEMdata, ResultsDataKeeper &RESdata, std::string output_vtk, bool PRINT_DEBUG_INFO) {
  /*
  std::string path_stress = FEMdata.res_dir + "/output/out_stress_" + FEMdata.get_name() + ".txt";
  if (PRINT_DEBUG_INFO) {
    std::cout << "StressComponents Size = " << RESdata.StressComponents.size() << "\n";
  }
  fstream out1;

  out1.open(path_stress, fstream::out);
  out1 << FEMdata.nodesCount << " " << FEMdata.elementsCount << "\n";
  for (int i = 0; i < RESdata.StressComponents.size(); i+=2) {
    out1 << RESdata.StressComponents[i] << " " << RESdata.StressComponents[i + 1] << "\n";
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
  */

  MakeVTKfile3D(output_vtk, FEMdata.nodesX, FEMdata.nodesY, FEMdata.nodesZ, FEMdata.elements,
              FEMdata.displacements, RESdata.Stress, RESdata.sigma_mises, RESdata.Deformation, RESdata.epsilon_mises);
}

//SparseMatrixCOO AssemblyStiffnessMatrix(FEMdataKeeper &FEMdata) {
//    vecSparseMatrixCOO triplets;
//    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
//        for (int i = 0; i < 4; ++i) {
//            for (int j = 0; j < 4; ++j) {
//                for (int ilocal = 0; ilocal < 3; ++ilocal) {
//                    for (int jlocal = 0; jlocal < 3; ++jlocal) {
//                        float value = it->Klocal(3 * i + ilocal, 3 * j + jlocal);
//                        if (value != 0.0) {
//                            Triplet tmp(3 * it->nodesIds[i] + ilocal, 3 * it->nodesIds[j] + jlocal, value);
//                            triplets.push_back(tmp);
//                        }
//                    }
//                }
//            }
//        }
//    }

//    cout << "CalculateStiffnessMatrix success\n";
//    cout << "Triplets Size = " << triplets.size() << std::endl;

//    SparseMatrixCOO globalK(triplets.size());
//    globalK.ConvertTripletToSparse(triplets);

//    std::cout << "new size= "<<globalK.get_size()<<"\n";

//    return globalK;
//}

void ApplyConstraints_EbE(FEMdataKeeper &FEMdata) {
  CheckRunTime(__func__)
  std::vector<int> indicesToConstraint;
  for (std::vector<Constraint>::const_iterator it = FEMdata.constraints.begin(); it != FEMdata.constraints.end(); ++it) {
    if (it->type & Constraint::UX) {
      indicesToConstraint.push_back(3 * it->node + 0);
    }
    if (it->type & Constraint::UY) {
      indicesToConstraint.push_back(3 * it->node + 1);
    }
    if (it->type & Constraint::UZ) {
      indicesToConstraint.push_back(3 * it->node + 2);
    }
  }

  //CUDA
  std::cout << "CONSTRAINTS\n";
  FEMdata.CudaIndicesToConstraints.Resize(indicesToConstraint.size());
  for (int i = 0; i < indicesToConstraint.size(); ++i) {
    FEMdata.CudaIndicesToConstraints[i] = indicesToConstraint[i];
//    std::cout << FEMdata.CudaIndicesToConstraints[i] << " ";
  }
  FEMdata.CudaIndicesToConstraintsCount = indicesToConstraint.size();
}

void gpuPCG_EbE_vec(FEMdataKeeper &FEMdata, MyArray &res, bool doAssemblyRes, float eps, bool PRINT_DEBUG_INFO) {
  CheckRunTime(__func__)
  int n_elems  = FEMdata.elementsCount;
  int n_gl_dofs = FEMdata.nodesCount * DIM;
  int grid_size = 6 * (DIM - 1) * n_elems;

  gpuDataKeeper gpu_data(FEMdata.elementsCount, FEMdata.nodesCount, doAssemblyRes);

  copyElementsAndFlocals(FEMdata, gpu_data);
  gpuGenerateMask(gpu_data, FEMdata.elementsCount);
  gpuCountNAdjElem(gpu_data, grid_size);

//  gpuCalculateKlocal2(gpu_data, FEMdata.elementsCount, FEMdata.nodesX.get_data(), FEMdata.nodesY.get_data(),
//                      FEMdata.nodesCount, FEMdata.D.get_data(), FEMdata.CudaIndicesToConstraints.get_data(),
//                      FEMdata.CudaIndicesToConstraintsCount);
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
    gpuMultiplyKlocalByVec(gpu_data, FEMdata.elementsCount);
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
//  res.Show();

}

void CalculateFEM_EbE_vec_GPU(FEMdataKeeper &FEMdata, bool PRINT_DEBUG_INFO) {
  CheckRunTime(__func__)
  //    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
  //        it->CalculateKlocal(FEMdata.D, FEMdata.nodesX, FEMdata.nodesY);
  //    }

  for (std::vector<BoundaryEdge>::iterator it = FEMdata.boundary.begin(); it != FEMdata.boundary.end(); ++it) {
    FEMdata.elements[it->adj_elem1].CalculateFlocal(*it, FEMdata.nodesX, FEMdata.nodesY, FEMdata.nodesZ, 0.0f);
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

  gpuPCG_EbE_vec(FEMdata, FEMdata.displacements, true, 1e-4f, PRINT_DEBUG_INFO);
}
