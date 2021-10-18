
#include "femstruct.h"

using namespace std;

// CalculateStiffnessMatrix will be deprecated!
void Element::CalculateStiffnessMatrix(Matrix& D, std::vector<Triplet>& triplets, MyArray& nodesX, MyArray& nodesY, MyArray& nodesZ) {
    MyArray x(4), y(4), z(4);
    x[0] = nodesX[nodesIds[0]]; x[1] = nodesX[nodesIds[1]]; x[2] = nodesX[nodesIds[2]], x[3] = nodesX[nodesIds[3]];
    y[0] = nodesY[nodesIds[0]]; y[1] = nodesY[nodesIds[1]]; y[2] = nodesY[nodesIds[2]], y[3] = nodesY[nodesIds[3]];
    z[0] = nodesZ[nodesIds[0]]; z[1] = nodesZ[nodesIds[1]]; z[2] = nodesZ[nodesIds[2]], z[3] = nodesZ[nodesIds[3]];

    Matrix C(4, 4);
    C(0, 0) = C(1, 0) = C(2, 0) = C(3, 0) = 1.0;
    C(0, 1) = x[0]; C(1, 1) = x[1]; C(2, 1) = x[2]; C(3, 1) = x[3];
    C(0, 2) = y[0]; C(1, 2) = y[1]; C(2, 2) = y[2]; C(3, 2) = y[3];
    C(0, 3) = z[0]; C(1, 3) = z[1]; C(2, 3) = z[2]; C(3, 3) = z[3];

    Matrix IC(4, 4);
    C.inverse(IC, 4, 0);

    for (int i = 0; i < 4; i++) {
        B(0, 3 * i + 0) = IC(1, i);
        B(0, 3 * i + 1) = 0.0;
        B(0, 3 * i + 2) = 0.0;

        B(1, 3 * i + 0) = 0.0;
        B(1, 3 * i + 1) = IC(2, i);
        B(1, 3 * i + 2) = 0.0;

        B(2, 3 * i + 0) = 0.0;
        B(2, 3 * i + 1) = 0.0;
        B(2, 3 * i + 2) = IC(3, i);

        B(3, 3 * i + 0) = IC(2, i);
        B(3, 3 * i + 1) = IC(1, i);
        B(3, 3 * i + 2) = 0.0;

        B(4, 3 * i + 0) = 0.0;
        B(4, 3 * i + 1) = IC(3, i);
        B(4, 3 * i + 2) = IC(2, i);

        B(5, 3 * i + 0) = IC(3, i);
        B(5, 3 * i + 1) = 0.0;
        B(5, 3 * i + 2) = IC(1, i);
    }

    //Matrix K(12, 12);
    Matrix temp1(12, 6);
    Matrix temp_B(6, 12);
    float determinant = C.det(4);

    //K = B.transpose() * D * B * std::abs(C.det()) / 2.0;

    temp_B = B;

    temp_B.transpose2();

    temp1 = temp_B.Product(D);

    Klocal = temp1.Product(B);

    Klocal.scale(std::abs(determinant) / 6.0);

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            Triplet trplt11(3 * nodesIds[i] + 0, 3 * nodesIds[j] + 0, Klocal(3 * i + 0, 3 * j + 0));
            Triplet trplt12(3 * nodesIds[i] + 0, 3 * nodesIds[j] + 1, Klocal(3 * i + 0, 3 * j + 1));
            Triplet trplt13(3 * nodesIds[i] + 0, 3 * nodesIds[j] + 2, Klocal(3 * i + 0, 3 * j + 2));

            Triplet trplt21(3 * nodesIds[i] + 1, 3 * nodesIds[j] + 0, Klocal(3 * i + 1, 3 * j + 0));
            Triplet trplt22(3 * nodesIds[i] + 1, 3 * nodesIds[j] + 1, Klocal(3 * i + 1, 3 * j + 1));
            Triplet trplt23(3 * nodesIds[i] + 1, 3 * nodesIds[j] + 2, Klocal(3 * i + 1, 3 * j + 2));

            Triplet trplt31(3 * nodesIds[i] + 2, 3 * nodesIds[j] + 0, Klocal(3 * i + 2, 3 * j + 0));
            Triplet trplt32(3 * nodesIds[i] + 2, 3 * nodesIds[j] + 1, Klocal(3 * i + 2, 3 * j + 1));
            Triplet trplt33(3 * nodesIds[i] + 2, 3 * nodesIds[j] + 2, Klocal(3 * i + 2, 3 * j + 2));

            if (trplt11.get_value() != 0.0) {
                triplets.push_back(trplt11);
            }
            if (trplt12.get_value() != 0.0) {
                triplets.push_back(trplt12);
            }
            if (trplt13.get_value() != 0.0) {
                triplets.push_back(trplt13);
            }
            if (trplt21.get_value() != 0.0) {
                triplets.push_back(trplt21);
            }
            if (trplt22.get_value() != 0.0) {
                triplets.push_back(trplt22);
            }
            if (trplt23.get_value() != 0.0) {
                triplets.push_back(trplt23);
            }
            if (trplt31.get_value() != 0.0) {
                triplets.push_back(trplt31);
            }
            if (trplt32.get_value() != 0.0) {
                triplets.push_back(trplt32);
            }
            if (trplt33.get_value() != 0.0) {
                triplets.push_back(trplt33);
            }
        }
    }
}

void Element::CalculateKlocal(Matrix& D, MyArray& nodesX, MyArray& nodesY, MyArray& nodesZ)
{
    MyArray x(4), y(4), z(4);
    x[0] = nodesX[nodesIds[0]]; x[1] = nodesX[nodesIds[1]]; x[2] = nodesX[nodesIds[2]], x[3] = nodesX[nodesIds[3]];
    y[0] = nodesY[nodesIds[0]]; y[1] = nodesY[nodesIds[1]]; y[2] = nodesY[nodesIds[2]], y[3] = nodesY[nodesIds[3]];
    z[0] = nodesZ[nodesIds[0]]; z[1] = nodesZ[nodesIds[1]]; z[2] = nodesZ[nodesIds[2]], z[3] = nodesZ[nodesIds[3]];

    Matrix C(4, 4);
    C(0, 0) = C(1, 0) = C(2, 0) = C(3, 0) = 1.0;
    C(0, 1) = x[0]; C(1, 1) = x[1]; C(2, 1) = x[2]; C(3, 1) = x[3];
    C(0, 2) = y[0]; C(1, 2) = y[1]; C(2, 2) = y[2]; C(3, 2) = y[3];
    C(0, 3) = z[0]; C(1, 3) = z[1]; C(2, 3) = z[2]; C(3, 3) = z[3];

    Matrix IC(4, 4);
    C.inverse(IC, 4, 0);

    for (int i = 0; i < 4; i++) {
        B(0, 3 * i + 0) = IC(1, i);
        B(0, 3 * i + 1) = 0.0;
        B(0, 3 * i + 2) = 0.0;

        B(1, 3 * i + 0) = 0.0;
        B(1, 3 * i + 1) = IC(2, i);
        B(1, 3 * i + 2) = 0.0;

        B(2, 3 * i + 0) = 0.0;
        B(2, 3 * i + 1) = 0.0;
        B(2, 3 * i + 2) = IC(3, i);

        B(3, 3 * i + 0) = IC(2, i);
        B(3, 3 * i + 1) = IC(1, i);
        B(3, 3 * i + 2) = 0.0;

        B(4, 3 * i + 0) = 0.0;
        B(4, 3 * i + 1) = IC(3, i);
        B(4, 3 * i + 2) = IC(2, i);

        B(5, 3 * i + 0) = IC(3, i);
        B(5, 3 * i + 1) = 0.0;
        B(5, 3 * i + 2) = IC(1, i);
    }

    //Matrix K(12, 12);
    Matrix temp1(12, 6);
    Matrix temp_B(6, 12);
    float determinant = C.det(4);

    //K = B.transpose() * D * B * std::abs(C.det()) / 2.0;

    temp_B = B;

    temp_B.transpose2();

    temp1 = temp_B.Product(D);

    Klocal = temp1.Product(B);

    Klocal.scale(std::abs(determinant) / 6.0);
}
