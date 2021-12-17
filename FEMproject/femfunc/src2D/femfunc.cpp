
#include "femfunc.h"

using namespace std;


float SetConstraints(int i, int j, float v, int index) {
    if (i == index || j == index) {
        return i == j ? 1.0 : 0.0;
    } else {
        return v;
    }
}

void ApplyConstraints(SparseMatrixCOO& K, MyArray& F, const std::vector<Constraint>& constraints, int n) {
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

void CalculateStressAndDeformation(std::vector<MyArray> &Deformation,
                                   std::vector<MyArray> &Stress,
                                   std::vector<float> &epsilon_mises,
                                   std::vector<float> &sigma_mises,
                                   Matrix D,
                                   std::vector<Element> elements,
                                   MyArray displacements) {
    CheckRunTime(__func__)
    MyArray StressVector(3);
    MyArray DeformationVector(3);
    MyArray delta(6);

    for (std::vector<Element>::iterator it = elements.begin(); it != elements.end(); ++it) {
        delta[0] = displacements[2 * it->nodesIds[0] + 0];
        delta[1] = displacements[2 * it->nodesIds[0] + 1];
        delta[2] = displacements[2 * it->nodesIds[1] + 0];
        delta[3] = displacements[2 * it->nodesIds[1] + 1];
        delta[4] = displacements[2 * it->nodesIds[2] + 0];
        delta[5] = displacements[2 * it->nodesIds[2] + 1];

        DeformationVector = it->B.Product(delta);
        StressVector = D.Product(DeformationVector);

        double sigma = sqrt(StressVector[0] * StressVector[0] - StressVector[0]
                        * StressVector[1] + StressVector[1] * StressVector[1] + 3.0 * StressVector[2] * StressVector[2]);
        sigma_mises.push_back(sigma);

        double epsilon = sqrt(DeformationVector[0] * DeformationVector[0] - DeformationVector[0]
                        * DeformationVector[1] + DeformationVector[1] * DeformationVector[1] + 3.0 * DeformationVector[2] * DeformationVector[2]);
        epsilon_mises.push_back(epsilon);

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

    return false;
}

void CalculateStressAlongAxe(std::vector<float> &StressComponents,
                             std::string axe,
                             std::string stress_component,
                             float fixed_value,
                             float a,
                             float b,
                             std::vector<MyArray> Stress,
                             MyArray nodesX,
                             MyArray nodesY,
                             std::vector<Element> elements) {
    CheckRunTime(__func__)
    MyArray range(100);
    float h = (b - a) / (range.get_size() - 1);
    range[0] = a;
    for (int i = 1; i < range.get_size(); ++i) {
        range[i] = range[i - 1] + h;
    }
    int component_id = 0;
    if (stress_component == "x") {
        component_id = 0;
    } else if (stress_component == "y") {
        component_id = 1;
    } else if (stress_component == "xy") {
        component_id = 2;
    }

    for (int i = 0; i < range.get_size(); ++i) {
        float x = (axe == "x") ? range[i] : fixed_value;
        float y = (axe == "y") ? range[i] : fixed_value;
        for (int j = 0; j < elements.size(); ++j) {
            float x1 = nodesX[elements[j].nodesIds[0]], y1 = nodesY[elements[j].nodesIds[0]];
            float x2 = nodesX[elements[j].nodesIds[1]], y2 = nodesY[elements[j].nodesIds[1]];
            float x3 = nodesX[elements[j].nodesIds[2]], y3 = nodesY[elements[j].nodesIds[2]];
            if (CheckPointInside(x, y, x1, y1, x2, y2, x3, y3)) {
                StressComponents.push_back(range[i]);
                StressComponents.push_back(Stress[j][component_id]);
                break;
            }
        }
    }
}

void CalculateFiniteElementMethod(FEMdataKeeper &FEMdata) {
    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
        it->CalculateKlocal(FEMdata.D, FEMdata.nodesX, FEMdata.nodesY);
    }

    int num = 0;
    for (std::vector<BoundaryEdge>::iterator it = FEMdata.boundary.begin(); it != FEMdata.boundary.end(); ++it) {
        FEMdata.elements[it->adj_elem1].CalculateFlocal(*it, FEMdata.nodesX, FEMdata.nodesY, FEMdata.pressure[num++]);
    }

    SparseMatrixCOO globalK = AssemblyStiffnessMatrix   (FEMdata);
    //globalK.resize();

    //SparseMatrixCOO globalK2(nonzero);
    //globalK2 = globalK.DeleteZeros();

    MyArray F = AssemblyF(FEMdata); // globalK * displacements = F
    F.add(FEMdata.loads);

    ApplyConstraints(globalK, F, FEMdata.constraints, FEMdata.nodesCount);

    globalK.SortIt();

    int nonzero = globalK.CountNonZero();
    cout << "nonzero = " << nonzero << endl;

    //globalK.ShowAsMatrixNumber(0, globalK.get_size(), 2*FEMdata.nodesCount);

//    Matrix temp(2*FEMdata.nodesCount);
//    globalK.ConvertToMatrix(temp);
//    temp.LU_solve(temp, F, FEMdata.displacements, 2*FEMdata.nodesCount);

    globalK.CGM_solve(F, FEMdata.displacements, 1e-5);
}

void MakeResults(FEMdataKeeper &FEMdata, std::string output_vtk) {
    //POSTPROCESSING
    std::vector<MyArray> Deformation;
    std::vector<MyArray> Stress;
    std::vector<float> sigma_mises;
    std::vector<float> epsilon_mises;
    std::vector<float> StressComponents;

    CalculateStressAndDeformation(Deformation, Stress, epsilon_mises, sigma_mises, FEMdata.D, FEMdata.elements, FEMdata.displacements);

//    CalculateStressAlongAxe(StressComponents, "x", "xy", -2.8, -4.0, 2.0, Stress, FEMdata.nodesX, FEMdata.nodesY, FEMdata.elements);
//    std::string path = "C:/Users/mokin/Desktop/git/CudaFEMproject/prak_results/out_stress_ " + FEMdata.get_name() + ".txt";
//    std::cout << "StressComponents Size = " << StressComponents.size() << "\n";
//    fstream out;
//    out.open(path, fstream::out);
//    for (int i = 0; i < StressComponents.size(); i+=2) {
//        out << StressComponents[i] << " " << StressComponents[i + 1] << "\n";
//    }

    MakeVTKfile2D(output_vtk, FEMdata.nodesX, FEMdata.nodesY, FEMdata.elements,
                FEMdata.displacements, Stress, sigma_mises, Deformation, epsilon_mises);
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
    MyArray F(DIM * FEMdata.nodesCount);
    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
        for (int i = 0; i < 3; ++i) {
            F[2 * it->nodesIds[i] + 0] += it->Flocal[2 * i + 0];
            F[2 * it->nodesIds[i] + 1] += it->Flocal[2 * i + 1];
        }
    }
    return F;
}
