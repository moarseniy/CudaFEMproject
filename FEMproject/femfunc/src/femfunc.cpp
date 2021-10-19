
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
                                   Matrix D,
                                   std::vector<Element> elements,
                                   MyArray displacements) {
    CheckRunTime(__func__)
    MyArray StressVector(6);
    MyArray DeformationVector(6);
    MyArray delta(12);

    for (std::vector<Element>::iterator it = elements.begin(); it != elements.end(); ++it) {
        delta[0] = displacements[3 * it->nodesIds[0] + 0];
        delta[1] = displacements[3 * it->nodesIds[0] + 1];
        delta[2] = displacements[3 * it->nodesIds[0] + 2];
        delta[3] = displacements[3 * it->nodesIds[1] + 0];
        delta[4] = displacements[3 * it->nodesIds[1] + 1];
        delta[5] = displacements[3 * it->nodesIds[1] + 2];
        delta[6] = displacements[3 * it->nodesIds[2] + 0];
        delta[7] = displacements[3 * it->nodesIds[2] + 1];
        delta[8] = displacements[3 * it->nodesIds[2] + 2];
        delta[9] = displacements[3 * it->nodesIds[3] + 0];
        delta[10] = displacements[3 * it->nodesIds[3] + 1];
        delta[11] = displacements[3 * it->nodesIds[3] + 2];

        DeformationVector = it->B.Product(delta);
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

void CalculateFiniteElementMethod(FEMdataKeeper &FEMdata) {
    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
        it->CalculateKlocal(FEMdata.D, FEMdata.nodesX, FEMdata.nodesY, FEMdata.nodesZ);
    }

    SparseMatrixCOO globalK = AssemblyStiffnessMatrix(FEMdata);

    ApplyConstraints(globalK, FEMdata.constraints, FEMdata.loads.get_size());

    globalK.SortIt();

    cout << "nonzero = " << globalK.CountNonZero() << endl;

    globalK.CGM_solve(FEMdata.loads, FEMdata.displacements, FEMdata.loads.get_size(), 1e-10);
}

void MakeResults(FEMdataKeeper &FEMdata, std::string output_vtk) {
    //POSTPROCESSING
    std::vector<MyArray> Deformation;
    std::vector<MyArray> Stress;
    std::vector<float> sigma_mises;
    std::vector<float> epsilon_mises;

    CalculateStressAndDeformation(Deformation, Stress, epsilon_mises, sigma_mises, FEMdata.D, FEMdata.elements, FEMdata.displacements);

    MakeVTKfile(output_vtk, FEMdata.nodesX, FEMdata.nodesY, FEMdata.nodesZ, FEMdata.elements,
                FEMdata.displacements, Stress, sigma_mises, Deformation, epsilon_mises);
}

SparseMatrixCOO AssemblyStiffnessMatrix(FEMdataKeeper &FEMdata)
{
    vecSparseMatrixCOO triplets;
    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            for (int ilocal = 0; ilocal < 3; ++ilocal) {
                for (int jlocal = 0; jlocal < 3; ++jlocal) {
                    float value = it->Klocal(3 * i + ilocal, 3 * j + jlocal);
                    if (value != 0.0) {
                        Triplet tmp(3 * it->nodesIds[i] + ilocal, 3 * it->nodesIds[j] + jlocal, value);
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
