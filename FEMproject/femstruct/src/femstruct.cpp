
#include "femstruct.h"

using namespace std;

// CalculateStiffnessMatrix will be deprecated!
//void Element::CalculateStiffnessMatrix(Matrix& D, std::vector<Triplet>& triplets, MyArray& nodesX, MyArray& nodesY, MyArray& nodesZ) {
//    MyArray x(4), y(4), z(4);
//    x[0] = nodesX[nodesIds[0]]; x[1] = nodesX[nodesIds[1]]; x[2] = nodesX[nodesIds[2]], x[3] = nodesX[nodesIds[3]];
//    y[0] = nodesY[nodesIds[0]]; y[1] = nodesY[nodesIds[1]]; y[2] = nodesY[nodesIds[2]], y[3] = nodesY[nodesIds[3]];
//    z[0] = nodesZ[nodesIds[0]]; z[1] = nodesZ[nodesIds[1]]; z[2] = nodesZ[nodesIds[2]], z[3] = nodesZ[nodesIds[3]];

//    Matrix C(4, 4);
//    C(0, 0) = C(1, 0) = C(2, 0) = C(3, 0) = 1.0;
//    C(0, 1) = x[0]; C(1, 1) = x[1]; C(2, 1) = x[2]; C(3, 1) = x[3];
//    C(0, 2) = y[0]; C(1, 2) = y[1]; C(2, 2) = y[2]; C(3, 2) = y[3];
//    C(0, 3) = z[0]; C(1, 3) = z[1]; C(2, 3) = z[2]; C(3, 3) = z[3];

//    Matrix IC(4, 4);
//    C.inverse(IC, 4, 0);

//    for (int i = 0; i < 4; i++) {
//        B(0, 3 * i + 0) = IC(1, i);
//        B(0, 3 * i + 1) = 0.0;
//        B(0, 3 * i + 2) = 0.0;

//        B(1, 3 * i + 0) = 0.0;
//        B(1, 3 * i + 1) = IC(2, i);
//        B(1, 3 * i + 2) = 0.0;

//        B(2, 3 * i + 0) = 0.0;
//        B(2, 3 * i + 1) = 0.0;
//        B(2, 3 * i + 2) = IC(3, i);

//        B(3, 3 * i + 0) = IC(2, i);
//        B(3, 3 * i + 1) = IC(1, i);
//        B(3, 3 * i + 2) = 0.0;

//        B(4, 3 * i + 0) = 0.0;
//        B(4, 3 * i + 1) = IC(3, i);
//        B(4, 3 * i + 2) = IC(2, i);

//        B(5, 3 * i + 0) = IC(3, i);
//        B(5, 3 * i + 1) = 0.0;
//        B(5, 3 * i + 2) = IC(1, i);
//    }

//    //Matrix K(12, 12);
//    Matrix temp1(12, 6);
//    Matrix temp_B(6, 12);
//    float determinant = C.det(4);

//    //K = B.transpose() * D * B * std::abs(C.det()) / 2.0;

//    temp_B = B;

//    temp_B.transpose2();

//    temp1 = temp_B.Product(D);

//    Klocal = temp1.Product(B);

//    Klocal.scale(std::abs(determinant) / 6.0);

//    for (int i = 0; i < 4; i++) {
//        for (int j = 0; j < 4; j++) {
//            Triplet trplt11(3 * nodesIds[i] + 0, 3 * nodesIds[j] + 0, Klocal(3 * i + 0, 3 * j + 0));
//            Triplet trplt12(3 * nodesIds[i] + 0, 3 * nodesIds[j] + 1, Klocal(3 * i + 0, 3 * j + 1));
//            Triplet trplt13(3 * nodesIds[i] + 0, 3 * nodesIds[j] + 2, Klocal(3 * i + 0, 3 * j + 2));

//            Triplet trplt21(3 * nodesIds[i] + 1, 3 * nodesIds[j] + 0, Klocal(3 * i + 1, 3 * j + 0));
//            Triplet trplt22(3 * nodesIds[i] + 1, 3 * nodesIds[j] + 1, Klocal(3 * i + 1, 3 * j + 1));
//            Triplet trplt23(3 * nodesIds[i] + 1, 3 * nodesIds[j] + 2, Klocal(3 * i + 1, 3 * j + 2));

//            Triplet trplt31(3 * nodesIds[i] + 2, 3 * nodesIds[j] + 0, Klocal(3 * i + 2, 3 * j + 0));
//            Triplet trplt32(3 * nodesIds[i] + 2, 3 * nodesIds[j] + 1, Klocal(3 * i + 2, 3 * j + 1));
//            Triplet trplt33(3 * nodesIds[i] + 2, 3 * nodesIds[j] + 2, Klocal(3 * i + 2, 3 * j + 2));

//            if (trplt11.get_value() != 0.0) {
//                triplets.push_back(trplt11);
//            }
//            if (trplt12.get_value() != 0.0) {
//                triplets.push_back(trplt12);
//            }
//            if (trplt13.get_value() != 0.0) {
//                triplets.push_back(trplt13);
//            }
//            if (trplt21.get_value() != 0.0) {
//                triplets.push_back(trplt21);
//            }
//            if (trplt22.get_value() != 0.0) {
//                triplets.push_back(trplt22);
//            }
//            if (trplt23.get_value() != 0.0) {
//                triplets.push_back(trplt23);
//            }
//            if (trplt31.get_value() != 0.0) {
//                triplets.push_back(trplt31);
//            }
//            if (trplt32.get_value() != 0.0) {
//                triplets.push_back(trplt32);
//            }
//            if (trplt33.get_value() != 0.0) {
//                triplets.push_back(trplt33);
//            }
//        }
//    }
//}

void Element::CalculateKlocal(Matrix& D, MyArray& nodesX, MyArray& nodesY, MyArray& nodesZ) {
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

void Element::CalculateFlocal(BoundaryEdge& edge, MyArray& nodesX, MyArray& nodesY, MyArray& nodesZ, float t) {
  //CheckRunTime(__func__)
  edge.update(t);
  float pressure_value = edge.value;
  float X0 = nodesX[edge.node0], X1 = nodesX[edge.node1], X2 = nodesX[edge.node2];
  float Y0 = nodesY[edge.node0], Y1 = nodesY[edge.node1], Y2 = nodesY[edge.node2];
  float Z0 = nodesZ[edge.node0], Z1 = nodesZ[edge.node1], Z2 = nodesZ[edge.node2];
//  float edge_length = std::sqrt((X1 - X0) * (X1 - X0) + (Y1 - Y0) * (Y1 - Y0) + (Z1 - Z0) * (Z1 - Z0));
  struct MathVector{float x, y, z;};
  MathVector a{X1 - X0, Y1 - Y0, Z1 - Z0}, b{X2 - X0, Y2 - Y0, Z2 - Z0};
  float face_area = 0.5f * ((a.y * b.z - a.z * b.y) * (a.y * b.z - a.z * b.y) +
                            (a.z * b.x - a.x * b.z) * (a.z * b.x - a.x * b.z) +
                            (a.x * b.y - a.y * b.x) * (a.x * b.y - a.y * b.x));

  // Grisha: Consider adding the third node when parsing to avoid this abundance of if-else statements
  int edge2elem_num[4];
  if (edge.node0 == nodesIds[0]) {
    edge2elem_num[0] = 0;
  } else if (edge.node0 == nodesIds[1]) {
    edge2elem_num[0] = 1;
  } else if (edge.node0 == nodesIds[2]) {
    edge2elem_num[0] = 2;
  } else if (edge.node0 == nodesIds[3]) {
    edge2elem_num[0] = 3;
  }

  if (edge.node1 == nodesIds[0]) {
    edge2elem_num[1] = 0;
  } else if (edge.node1 == nodesIds[1]) {
    edge2elem_num[1] = 1;
  } else if (edge.node1 == nodesIds[2]) {
    edge2elem_num[1] = 2;
  } else if (edge.node1 == nodesIds[3]) {
    edge2elem_num[1] = 3;
  }

  if (edge.node2 == nodesIds[0]) {
    edge2elem_num[2] = 0;
  } else if (edge.node2 == nodesIds[1]) {
    edge2elem_num[2] = 1;
  } else if (edge.node2 == nodesIds[2]) {
    edge2elem_num[2] = 2;
  } else if (edge.node2 == nodesIds[3]) {
    edge2elem_num[2] = 3;
  }

  if (edge2elem_num[0] == 0) {              // 0___
    if (edge2elem_num[1] == 1) {            // 01__
        if (edge2elem_num[2] == 2) {        // 012_
            edge2elem_num[3] = 3;           // 0123
        } else {                            // 013_
            edge2elem_num[3] = 2;           // 0132
        }
    } else {
        if (edge2elem_num[1] == 2) {        // 02__
            if (edge2elem_num[2] == 1) {    // 021_
                edge2elem_num[2] = 3;       // 0213
            } else {                        // 023_
                edge2elem_num[2] = 1;       // 0231
            }
        } else {                            // 03__
            if (edge2elem_num[2] == 1) {    // 031_
                edge2elem_num[2] = 2;       // 0312
            } else {                        // 032_
                edge2elem_num[2] = 1;       // 0321
            }
        }
    }
  } else if (edge2elem_num[0] == 1) {       // 1___
    if (edge2elem_num[1] == 0) {            // 10__
        if (edge2elem_num[2] == 2) {        // 102_
            edge2elem_num[3] = 3;           // 1023
        } else {                            // 103_
            edge2elem_num[3] = 2;           // 1032
        }
    } else {
        if (edge2elem_num[1] == 2) {        // 12__
            if (edge2elem_num[2] == 0) {    // 120_
                edge2elem_num[2] = 3;       // 1203
            } else {                        // 123_
                edge2elem_num[2] = 0;       // 1230
            }
        } else {                            // 13__
            if (edge2elem_num[2] == 0) {    // 130_
                edge2elem_num[2] = 2;       // 1302
            } else {                        // 132_
                edge2elem_num[2] = 0;       // 1320
            }
        }
    }
  } else if (edge2elem_num[0] == 2) {       // 2___
    if (edge2elem_num[1] == 0) {            // 20__
        if (edge2elem_num[2] == 1) {        // 201_
            edge2elem_num[3] = 3;           // 2013
        } else {                            // 203_
            edge2elem_num[3] = 1;           // 2031
        }
    } else {
        if (edge2elem_num[1] == 1) {        // 21__
            if (edge2elem_num[2] == 0) {    // 210_
                edge2elem_num[2] = 3;       // 2103
            } else {                        // 213_
                edge2elem_num[2] = 0;       // 2130
            }
        } else {                            // 23__
            if (edge2elem_num[2] == 0) {    // 230_
                edge2elem_num[2] = 1;       // 2301
            } else {                        // 231_
                edge2elem_num[2] = 0;       // 2310
            }
        }
    }
  } else if (edge2elem_num[0] == 3) {       // 3___
    if (edge2elem_num[1] == 0) {            // 30__
        if (edge2elem_num[2] == 1) {        // 301_
            edge2elem_num[3] = 2;           // 3012
        } else {                            // 302_
            edge2elem_num[3] = 1;           // 3021
        }
    } else {
        if (edge2elem_num[1] == 1) {        // 31__
            if (edge2elem_num[2] == 0) {    // 310_
                edge2elem_num[2] = 2;       // 3102
            } else {                        // 312_
                edge2elem_num[2] = 0;       // 3120
            }
        } else {                            // 32__
            if (edge2elem_num[2] == 0) {    // 320_
                edge2elem_num[2] = 1;       // 3201
            } else {                        // 321_
                edge2elem_num[2] = 0;       // 3210
            }
        }
    }
  }

  if (edge.type0 & Constraint::UX) {
    Flocal[3 * edge2elem_num[0] + 0] = 0.0f;
  } else {
    Flocal[3 * edge2elem_num[0] + 0] = - 0.5 * pressure_value * face_area * edge.normal_x;
  }
  if (edge.type0 & Constraint::UY) {
    Flocal[3 * edge2elem_num[0] + 1] = 0.0f;
  } else {
    Flocal[3 * edge2elem_num[0] + 1] = - 0.5 * pressure_value * face_area * edge.normal_y;
  }
  if (edge.type0 & Constraint::UZ) {
    Flocal[3 * edge2elem_num[0] + 2] = 0.0f;
  } else {
    Flocal[3 * edge2elem_num[0] + 2] = - 0.5 * pressure_value * face_area * edge.normal_z;
  }

  if (edge.type1 & Constraint::UX) {
    Flocal[3 * edge2elem_num[1] + 0] = 0.0f;
  } else {
    Flocal[3 * edge2elem_num[1] + 0] = - 0.5 * pressure_value * face_area * edge.normal_x;
  }
  if (edge.type1 & Constraint::UY) {
    Flocal[3 * edge2elem_num[1] + 1] = 0.0f;
  } else {
    Flocal[3 * edge2elem_num[1] + 1] = - 0.5 * pressure_value * face_area * edge.normal_y;
  }
  if (edge.type1 & Constraint::UZ) {
    Flocal[3 * edge2elem_num[1] + 2] = 0.0f;
  } else {
    Flocal[3 * edge2elem_num[1] + 2] = - 0.5 * pressure_value * face_area * edge.normal_z;
  }

  if (edge.type2 & Constraint::UX) {
    Flocal[3 * edge2elem_num[2] + 0] = 0.0f;
  } else {
    Flocal[3 * edge2elem_num[2] + 0] = - 0.5 * pressure_value * face_area * edge.normal_x;
  }
  if (edge.type2 & Constraint::UY) {
    Flocal[3 * edge2elem_num[2] + 1] = 0.0f;
  } else {
    Flocal[3 * edge2elem_num[2] + 1] = - 0.5 * pressure_value * face_area * edge.normal_y;
  }
  if (edge.type2 & Constraint::UZ) {
    Flocal[3 * edge2elem_num[2] + 2] = 0.0f;
  } else {
    Flocal[3 * edge2elem_num[2] + 2] = - 0.5 * pressure_value * face_area * edge.normal_z;
  }

}

void Load::assignElement(std::unordered_map <int, std::vector<int>> nodeAdjElem) {
  this->elem = nodeAdjElem[this->dof / DIM][0];
}

void TimeDependentEntity::update(float t) {
  if (t < this->timeshift) {
    this->value = 0.0f;
  } else {
    (this->*wavelet)(t - this->timeshift);
  }
}

void TimeDependentEntity::Constant(float t) {
  this->value = ampl;
}

void TimeDependentEntity::Ricker(float t) {
  float tmp = M_PI * this->freq * t - M_PI;
  this->value = (1.0f - 2.0f * tmp*tmp) *
      std::exp(-1.0f * tmp*tmp);

  this->value *= this->ampl;
}

// Aldridge, D. F. (1990). The Berlage wavelet. GEOPHYSICS, 55(11), 1508–1511. doi:10.1190/1.1442799
// CAE-Fidesys-4.0/preprocessor/bin/help/finite_element_model/non_exodus/time_formulas.htm
void TimeDependentEntity::Berlage(float t) {
  float w0 = 2 * M_PI * this->freq;
  float w1 = w0 / sqrtf(3.0f);
  this->value = w1*w1/4.0f*std::exp(-1.0f*w1*t)*(
        std::sin(w0*t) * (1.0f/(w1*w1*w1) + t/(w1*w1) - t*t/w1) -
        std::cos(w0*t) * sqrtf(3.0f) * (t*t/w1 + t/(w1*w1)) );

  this->value *= this->ampl;
}

// I always liked strtof since it lets you specify an end pointer. [ https://stackoverflow.com/a/57163016 ]
bool TimeDependentEntity::isFloat(const std::string& str) {
  char* ptr;
  strtof(str.c_str(), &ptr);
  return (*ptr) == '\0';
}

void TimeDependentEntity::parseString(std::string& str) {
  if (this->isFloat(str)) {
    float str_value = std::stof(str);
    this->value = str_value;
    this->ampl = str_value;
  } else {             // if str is a time-dependent function, not constant value
    // parse string
    // https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c

    std::string wavelet_str = str.substr(0, str.find("("));
    str.erase(0, str.find("(") + 1);
    //std::cout << "DEBUG : xstr = " << str << std::endl;
    //std::cout << "DEBUG : wavelet_str = " << wavelet_str << std::endl;

    std::string freq_str = str.substr(0, str.find(","));
    str.erase(0, str.find(",") + 1);
    //std::cout << "DEBUG : xstr = " << str << std::endl;
    //std::cout << "DEBUG : freq_str = " << freq_str << std::endl;

    std::string timeshift_str = str.substr(0, str.find(","));
    str.erase(0, str.find(",") + 1);
    //std::cout << "DEBUG : xstr = " << str << std::endl;
    //std::cout << "DEBUG : timeshift_str = " << timeshift_str << std::endl;

    std::string ampl_str = str.substr(0, str.find(")"));
    //xstr.erase(0, xstr.find(")") + 1);
    //std::cout << "DEBUG : xstr = " << str << std::endl;
    //std::cout << "DEBUG : ampl_str = " << ampl_str << std::endl;

    this->value = std::stof(ampl_str);
    this->ampl  = std::stof(ampl_str);

    // Try reading if you want more elegancy:
    // https://stackoverflow.com/questions/650162/why-the-switch-statement-cannot-be-applied-on-strings

    // https://stackoverflow.com/a/313990
    std::transform(wavelet_str.begin(), wavelet_str.end(), wavelet_str.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    if (wavelet_str == "ricker") {
      this->wavelet = &TimeDependentEntity::Ricker;
    } else if (wavelet_str == "berlage") {
      this->wavelet = &TimeDependentEntity::Berlage;
    }

    //std::cout << "DEBUG : wavelet_str = " << wavelet_str << std::endl;

    this->timeshift = std::stof(timeshift_str);
    this->freq      = std::stof(freq_str);
    //std::cout << "value = " << this->value << "; timeshift = " << this->timeshift << "; load.freq = " << this->freq << std::endl;

  }
}
