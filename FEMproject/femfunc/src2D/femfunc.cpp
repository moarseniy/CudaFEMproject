
#include "femfunc.h"

using namespace std;


float SetConstraints(int i, int j, float v, int index) {
    if (i == index || j == index) {
        return i == j ? 1.0f : 0.0f;
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

void CalculateStressAlongAxis(std::vector<float> &StressComponents,
                             std::string axe,
                             std::string stress_component,
                             float fixed_value,
                             float a,
                             float b,
                             std::vector<MyArray> &Stress,
                             MyArray nodesX,
                             MyArray nodesY,
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
            float x1 = nodesX[elements[j].nodesIds[0]], y1 = nodesY[elements[j].nodesIds[0]];
            float x2 = nodesX[elements[j].nodesIds[1]], y2 = nodesY[elements[j].nodesIds[1]];
            float x3 = nodesX[elements[j].nodesIds[2]], y3 = nodesY[elements[j].nodesIds[2]];
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

void CalculateStressAlongAxisSmooth(std::vector<float> &StressComponentsSmooth,
                             std::string axe,
                             float fixed_value,
                             float a,
                             float b,
                             MyArray StressSmooth,
                             MyArray nodesX,
                             MyArray nodesY,
                             std::vector<Element> elements) {
    CheckRunTime(__func__)

    MyArray range(100, a, b);
    for (int i = 0; i < range.get_size(); ++i) {
        float x = (axe == "x") ? range[i] : fixed_value;
        float y = (axe == "y") ? range[i] : fixed_value;
        for (int j = 0; j < elements.size(); ++j) {
            float x1 = nodesX[elements[j].nodesIds[0]], y1 = nodesY[elements[j].nodesIds[0]];
            float x2 = nodesX[elements[j].nodesIds[1]], y2 = nodesY[elements[j].nodesIds[1]];
            float x3 = nodesX[elements[j].nodesIds[2]], y3 = nodesY[elements[j].nodesIds[2]];
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
                             MyArray nodesX,
                             MyArray nodesY,
                             std::vector<Element> elements) {
    CheckRunTime(__func__)

    MyArray range(100, a, b);
    for (int i = 0; i < range.get_size(); ++i) {
        float x = range[i];
        float y = k * x + m;
        for (int j = 0; j < elements.size(); ++j) {
            float x1 = nodesX[elements[j].nodesIds[0]], y1 = nodesY[elements[j].nodesIds[0]];
            float x2 = nodesX[elements[j].nodesIds[1]], y2 = nodesY[elements[j].nodesIds[1]];
            float x3 = nodesX[elements[j].nodesIds[2]], y3 = nodesY[elements[j].nodesIds[2]];
            if (CheckPointInside(x, y, x1, y1, x2, y2, x3, y3)) {
                MisesComponents.push_back(range[i]);
                MisesComponents.push_back(sigma_mises[j]);
                break;
            }
        }
    }
}

void CalculateNodeAdjElem(FEMdataKeeper FEMdata, std::unordered_map <int, std::vector<int>> &a) {
    int n = 0;
    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
        a[it->nodesIds[0]].push_back(n);
        a[it->nodesIds[1]].push_back(n);
        a[it->nodesIds[2]].push_back(n);
        n++;
    }
}

void ApplyConstraints2(FEMdataKeeper &FEMdata) {
    std::vector<int> indicesToConstraint;
    for (std::vector<Constraint>::const_iterator it = FEMdata.constraints.begin(); it != FEMdata.constraints.end(); ++it) {
        if (it->type & Constraint::UX) {
            indicesToConstraint.push_back(2 * it->node + 0);
        }
        if (it->type & Constraint::UY) {
            indicesToConstraint.push_back(2 * it->node + 1);
        }
    }

    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
        for (int c_id = 0; c_id < indicesToConstraint.size(); c_id++) {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    for (int ilocal = 0; ilocal < 2; ++ilocal) {
                        for (int jlocal = 0; jlocal < 2; ++jlocal) {
                            if (2 * it->nodesIds[i] + ilocal == indicesToConstraint[c_id] ||
                                2 * it->nodesIds[j] + jlocal == indicesToConstraint[c_id]) {
                                if (2 * it->nodesIds[i] + ilocal != 2 * it->nodesIds[j] + jlocal) {
                                    it->Klocal(2 * i + ilocal, 2 * j + jlocal) = 0.0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void ApplyConstraints_EbE(FEMdataKeeper &FEMdata, std::unordered_map <int, std::vector<int>> &nodeAdjElem) {
    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
        for (std::vector<Constraint>::const_iterator itc = FEMdata.constraints.begin(); itc != FEMdata.constraints.end(); ++itc) {
            for (int local_node = 0; local_node < 3; ++local_node) {
                if (it->nodesIds[local_node] == itc->node) {
                    std::vector<int> adjElems = nodeAdjElem[itc->node];
                    for (int i = 0; i < adjElems.size(); ++i) {
                        Element adj_elem = FEMdata.elements[adjElems[i]];

                        for (int adj_local_node = 0; adj_local_node < 3; ++adj_local_node) {
                            if (adj_elem.nodesIds[adj_local_node] == itc->node) {
//                                std::cout << "Before aplying constraints:" << std::endl;
//                                adj_elem.Klocal.Show();
                                if (itc->type & Constraint::UX) {
                                    int mIdx = 2 * adj_local_node + 0;

                                    adj_elem.Flocal[mIdx] = 0.0;

                                    adj_elem.Klocal(mIdx, (mIdx + 1) % 6) = 0.0;
                                    adj_elem.Klocal(mIdx, (mIdx + 2) % 6) = 0.0;
                                    adj_elem.Klocal(mIdx, (mIdx + 3) % 6) = 0.0;
                                    adj_elem.Klocal(mIdx, (mIdx + 4) % 6) = 0.0;
                                    adj_elem.Klocal(mIdx, (mIdx + 5) % 6) = 0.0;
                                    adj_elem.Klocal((mIdx + 1) % 6, mIdx) = 0.0;
                                    adj_elem.Klocal((mIdx + 2) % 6, mIdx) = 0.0;
                                    adj_elem.Klocal((mIdx + 3) % 6, mIdx) = 0.0;
                                    adj_elem.Klocal((mIdx + 4) % 6, mIdx) = 0.0;
                                    adj_elem.Klocal((mIdx + 5) % 6, mIdx) = 0.0;
                                }
                                if (itc->type & Constraint::UY) {
                                    int mIdx = 2 * adj_local_node + 1;

                                    adj_elem.Flocal[mIdx] = 0.0;

                                    adj_elem.Klocal(mIdx, (mIdx + 1) % 6) = 0.0;
                                    adj_elem.Klocal(mIdx, (mIdx + 2) % 6) = 0.0;
                                    adj_elem.Klocal(mIdx, (mIdx + 3) % 6) = 0.0;
                                    adj_elem.Klocal(mIdx, (mIdx + 4) % 6) = 0.0;
                                    adj_elem.Klocal(mIdx, (mIdx + 5) % 6) = 0.0;
                                    adj_elem.Klocal((mIdx + 1) % 6, mIdx) = 0.0;
                                    adj_elem.Klocal((mIdx + 2) % 6, mIdx) = 0.0;
                                    adj_elem.Klocal((mIdx + 3) % 6, mIdx) = 0.0;
                                    adj_elem.Klocal((mIdx + 4) % 6, mIdx) = 0.0;
                                    adj_elem.Klocal((mIdx + 5) % 6, mIdx) = 0.0;
                                }

//                                std::cout << "After aplying constraints:" << std::endl;
//                                adj_elem.Klocal.Show();
//                                std::cout << std::endl << std::endl;

                                break;
                            }
                        }
                    }
                }
            }
        }
    }
}

void PCG_EbE(FEMdataKeeper &FEMdata, std::unordered_map <int, std::vector<int>> &node2adj_elem, float eps) {
    // see article "A distributed memory parallel element-by-element scheme based on Jacobi-conditioned conjugate gradient for 3D finite element analysis"
    // by Yaoru Liu, Weiyuan Zhou, Qiang Yang

    float eps_div = 1e-30f;
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

    std::cout.precision(16);
    std::cout << "gamma0\t\t\t= " << gamma0 << std::endl << std::endl;

    int n_iter = 0;
    do {
    ++n_iter;
    int enum_it = 0;
    std::cout << "Iteration #" << n_iter << std::endl;
    float sumElem = 0.0f;
    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
        // (1a-d)
        it->u = it->Klocal.Product(it->p);
        float tmp = it->p.dot_product(it->u);
        sumElem += tmp;
    }

    float alpha = gamma / sumElem;

    std::cout << "alpha (gamma / sumElem)\t= " << alpha << std::endl;
    std::cout << "alpha numerator (gamma)\t= " << gamma << std::endl;
    std::cout << "alpha denominator\t= " << sumElem << std::endl;

    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
        // (2a)
        it->x.add_weighted(it->p, 1.0f, alpha);
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

    std::cout << "beta\t\t\t= " << gamma_new / gamma << std::endl;
    std::cout << "gamma_new\t\t= " << gamma_new << std::endl;

    gamma = gamma_new;

    std::cout << std::endl;

    } while (1);
}

void CalculateFEM(FEMdataKeeper &FEMdata) {
    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
        it->CalculateKlocal(FEMdata.D, FEMdata.nodesX, FEMdata.nodesY);
    }

    int num = 0;
    for (std::vector<BoundaryEdge>::iterator it = FEMdata.boundary.begin(); it != FEMdata.boundary.end(); ++it) {
        FEMdata.elements[it->adj_elem1].CalculateFlocal(*it, FEMdata.nodesX, FEMdata.nodesY, FEMdata.pressure[num++]);
    }

    SparseMatrixCOO globalK = AssemblyStiffnessMatrix   (FEMdata);

    MyArray F = AssemblyF(FEMdata); // globalK * displacements = F
    F.add(FEMdata.loads);

    ApplyConstraints(globalK, F, FEMdata.constraints, FEMdata.nodesCount);

    globalK.SortIt();
    globalK.set_diag_elements();
    std::vector<float> diag_elem = globalK.get_diag_elements();
    int nonzero = globalK.CountNonZero();
    cout << "nonzero = " << nonzero << endl;

    globalK.PCG_solve(F, FEMdata.displacements, 1e-10f);
}

void CalculateFEM_EbE(FEMdataKeeper &FEMdata) {
    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
        it->CalculateKlocal(FEMdata.D, FEMdata.nodesX, FEMdata.nodesY);
    }

    int num = 0;
    for (std::vector<BoundaryEdge>::iterator it = FEMdata.boundary.begin(); it != FEMdata.boundary.end(); ++it) {
        FEMdata.elements[it->adj_elem1].CalculateFlocal(*it, FEMdata.nodesX, FEMdata.nodesY, FEMdata.pressure[num++]);
    }

    std::unordered_map <int, std::vector<int>> nodeAdjElem;
    CalculateNodeAdjElem(FEMdata, nodeAdjElem);

    ApplyConstraints2(FEMdata);
    PCG_EbE(FEMdata, nodeAdjElem, 1e-10f);
    AssemblyX(FEMdata, nodeAdjElem);
}

void SmoothResults(std::string stress_component, MyArray &SmoothStress, std::vector<MyArray> Stress,
                   int nodesCount, MyArray nodesX, MyArray nodesY, std::vector<Element> elements) {
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
        x1 = nodesX[elements[i].nodesIds[0]]; y1 = nodesY[elements[i].nodesIds[0]];
        x2 = nodesX[elements[i].nodesIds[1]]; y2 = nodesY[elements[i].nodesIds[1]];
        x3 = nodesX[elements[i].nodesIds[2]]; y3 = nodesY[elements[i].nodesIds[2]];

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

    CalculateStressAndDeformation(RESdata.Deformation,
                                  RESdata.Stress,
                                  RESdata.epsilon_mises,
                                  RESdata.sigma_mises,
                                  FEMdata.D, FEMdata.elements, FEMdata.displacements);

    float fixed_value = -3.8f;// x -3.0; y -3.0
    float a = -3.0f;// x -3.0 y -4.0
    float b = 5.0f;// x 5.0 y -3.0;

    CalculateStressAlongAxis(RESdata.StressComponents,
                             "x", "xx", fixed_value, a, b,
                             RESdata.Stress, FEMdata.nodesX, FEMdata.nodesY, FEMdata.elements);

    if (RESdata.withSmooth) {
        SmoothResults("xx", RESdata.SmoothStress, RESdata.Stress, FEMdata.nodesCount, FEMdata.nodesX, FEMdata.nodesY, FEMdata.elements);
        CalculateStressAlongAxisSmooth(RESdata.StressComponentsSmooth,
                             "x", fixed_value, a, b, RESdata.SmoothStress, FEMdata.nodesX, FEMdata.nodesY, FEMdata.elements);
    }

    a = -3.0f;
    b = 3.0f;
    float k = -0.6655f;
    float m = -0.0035f;

    if (RESdata.withMises) {
        CalculateMisesAlongLineMises(RESdata.MisesComponents,
                               k, m, a, b, RESdata.sigma_mises, FEMdata.nodesX, FEMdata.nodesY, FEMdata.elements);
    }
}

void WriteResults(FEMdataKeeper FEMdata, ResultsDataKeeper RESdata, std::string output_vtk) {
    std::string path_stress = FEMdata.res_dir + "output/out_stress_" + FEMdata.get_name() + ".txt";
    std::string path_stress_smooth = FEMdata.res_dir + "output/out_stress_" + FEMdata.get_name() + "_smooth.txt";
    std::string path_stress_mises = FEMdata.res_dir + "output/out_stress_" + FEMdata.get_name() + "_mises.txt";

    std::cout << "StressComponents Size = " << RESdata.StressComponents.size() << "\n";
    std::cout << "StressComponentsSmooth Size = " << RESdata.StressComponentsSmooth.size() << "\n";
    std::cout << "MisesComponents Size = " << RESdata.MisesComponents.size() << "\n";

    fstream out1, out2, out3;

    out1.open(path_stress, fstream::out);
    out1 << FEMdata.nodesCount << " " << FEMdata.elementsCount << "\n";
    for (int i = 0; i < RESdata.StressComponents.size(); i+=2) {
        out1 << RESdata.StressComponents[i] << " " << RESdata.StressComponents[i + 1] << "\n";
    }

    if (RESdata.withSmooth) {
        out2.open(path_stress_smooth, fstream::out);
        out2 << FEMdata.nodesCount << " " << FEMdata.elementsCount << "\n";
        for (int i = 0; i < RESdata.StressComponentsSmooth.size(); i+=2) {
            out2 << RESdata.StressComponentsSmooth[i] << " " << RESdata.StressComponentsSmooth[i + 1] << "\n";
        }
    }

    if (RESdata.withMises) {
        out3.open(path_stress_mises, fstream::out);
        out3 << FEMdata.nodesCount << " " << FEMdata.elementsCount << "\n";
        for (int i = 0; i < RESdata.MisesComponents.size(); i+=2) {
            out3 << RESdata.MisesComponents[i] << " " << RESdata.MisesComponents[i + 1] << "\n";
        }
    }

    MakeVTKfile2D(output_vtk, FEMdata.nodesX, FEMdata.nodesY, FEMdata.elements,
                FEMdata.displacements, RESdata.Stress, RESdata.sigma_mises, RESdata.Deformation, RESdata.epsilon_mises, RESdata.SmoothStress);
}
