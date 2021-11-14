
#include "femstruct.h"

using namespace std;

// CalculateStiffnessMatrix will be deprecated!
void Element::CalculateStiffnessMatrix(Matrix& D, std::vector<Triplet>& triplets, MyArray& nodesX, MyArray& nodesY) {
    MyArray x(3), y(3);
    x[0] = nodesX[nodesIds[0]]; x[1] = nodesX[nodesIds[1]]; x[2] = nodesX[nodesIds[2]];
    y[0] = nodesY[nodesIds[0]]; y[1] = nodesY[nodesIds[1]]; y[2] = nodesY[nodesIds[2]];

//    x.Show();
//    y.Show();
//    cout << endl;

    //Matrix B(3, 6);

    Matrix C(3, 3);
    C(0, 0) = C(1, 0) = C(2, 0) = 1.0;
    C(0, 1) = x[0]; C(1, 1) = x[1]; C(2, 1) = x[2];
    C(0, 2) = y[0]; C(1, 2) = y[1]; C(2, 2) = y[2];

    Matrix IC(3, 3);
    //C.Show();
    C.inverse(IC, 3, 0);

    for (int i = 0; i < 3; i++) {
        B(0, 2 * i + 0) = IC(1, i);
        B(0, 2 * i + 1) = 0.0;
        B(1, 2 * i + 0) = 0.0;
        B(1, 2 * i + 1) = IC(2, i);
        B(2, 2 * i + 0) = IC(2, i);
        B(2, 2 * i + 1) = IC(1, i);
    }
    //B.Show();

    Matrix K(6, 6);
    Matrix temp1(6, 3);
    Matrix temp2(6, 3);
    Matrix temp_B(3, 6);
    double determinant = C.det(3);

    //K = B.transpose() * D * B * std::abs(C.det()) / 2.0;

    temp_B = B;
    //temp_B.Show();

    temp_B.transpose2();
    //temp_B.Show();

    temp1 = temp_B.Product(D);
    //temp2.Show();

    K = temp1.Product(B);

    //K.Show();

    K.scale(std::abs(determinant) * 0.5);
    //K.Show();

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Triplet trplt11(2 * nodesIds[i] + 0, 2 * nodesIds[j] + 0, K(2 * i + 0, 2 * j + 0));
            Triplet trplt12(2 * nodesIds[i] + 0, 2 * nodesIds[j] + 1, K(2 * i + 0, 2 * j + 1));
            Triplet trplt21(2 * nodesIds[i] + 1, 2 * nodesIds[j] + 0, K(2 * i + 1, 2 * j + 0));
            Triplet trplt22(2 * nodesIds[i] + 1, 2 * nodesIds[j] + 1, K(2 * i + 1, 2 * j + 1));

            if (trplt11.get_value() != 0) {
                triplets.push_back(trplt11);
            }
            if (trplt12.get_value() != 0) {
                triplets.push_back(trplt12);
            }
            if (trplt21.get_value() != 0) {
                triplets.push_back(trplt21);
            }
            if (trplt22.get_value() != 0) {
                triplets.push_back(trplt22);
            }
        }
    }
}

void Element::CalculateKlocal(Matrix& D, MyArray& nodesX, MyArray& nodesY) {
    MyArray x(3), y(3);
    x[0] = nodesX[nodesIds[0]]; x[1] = nodesX[nodesIds[1]]; x[2] = nodesX[nodesIds[2]];
    y[0] = nodesY[nodesIds[0]]; y[1] = nodesY[nodesIds[1]]; y[2] = nodesY[nodesIds[2]];

//    x.Show();
//    y.Show();
//    cout << endl;

    //Matrix B(3, 6);

    Matrix C(3, 3);
    C(0, 0) = C(1, 0) = C(2, 0) = 1.0;
    C(0, 1) = x[0]; C(1, 1) = x[1]; C(2, 1) = x[2];
    C(0, 2) = y[0]; C(1, 2) = y[1]; C(2, 2) = y[2];

    Matrix IC(3, 3);
    //C.Show();
    C.inverse(IC, 3, 0);

    for (int i = 0; i < 3; i++) {
        B(0, 2 * i + 0) = IC(1, i);
        B(0, 2 * i + 1) = 0.0;
        B(1, 2 * i + 0) = 0.0;
        B(1, 2 * i + 1) = IC(2, i);
        B(2, 2 * i + 0) = IC(2, i);
        B(2, 2 * i + 1) = IC(1, i);
    }
    //B.Show();

    Matrix temp1(6, 3);
    Matrix temp2(6, 3);
    Matrix temp_B(3, 6);
    double determinant = C.det(3);

    //K = B.transpose() * D * B * std::abs(C.det()) / 2.0;

    temp_B = B;

    temp_B.transpose2();

    temp1 = temp_B.Product(D);

    Klocal = temp1.Product(B);

    Klocal.scale(std::abs(determinant) * 0.5);
}

