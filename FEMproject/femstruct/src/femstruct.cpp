#include <femstruct.h>

Element::Element(int dim) :
  DIM(dim),
  nodesIds(DIM + 1),
  B(3 * (DIM - 1), 6 * (DIM - 1)),
  Klocal(6 * (DIM - 1), 6 * (DIM - 1)),
  Mlocal(6 * (DIM - 1), 6 * (DIM - 1)),
  Clocal(6 * (DIM - 1), 6 * (DIM - 1)),
  Flocal(6 * (DIM - 1), 1),
  Alocal(6 * (DIM - 1), 6 * (DIM - 1)),
  blocal(6 * (DIM - 1), 1),
  res(6 * (DIM - 1), 1),
  x_pred(3 * DIM, 1),
  vel_pred(3 * DIM, 1),
  vel(3 * DIM, 1),
  m(3 * DIM, 1),
  x(3 * DIM, 1),
  r(3 * DIM, 1),
  z(3 * DIM, 1),
  s(3 * DIM, 1),
  p(3 * DIM, 1),
  u(3 * DIM, 1),
  b(3 * DIM, 1) {}

Element::~Element() {}

float Element::calculateArea(CPU_Matrix &Coordinates) {
  CPU_Matrix C(DIM + 1, DIM + 1);
  for (size_t i = 0; i < C.get_numRows(); ++i) {
    for (size_t j = 0; j < C.get_numCols(); ++j) {
      C.get_data()[j + i * C.get_numCols()] =
          j == 0 ? 1.f : Coordinates(j - 1, i); // 2 x 3
    }
  }
  return C.det();
}

void Element::calculateGradientMatrix(CPU_Matrix &Coordinates, float area) {
  B(0, 0) = Coordinates(1, 1) - Coordinates(1, 2);
  B(0, 1) = 0.f;
  B(0, 2) = Coordinates(1, 2) - Coordinates(1, 0);
  B(0, 3) = 0.f;
  B(0, 4) = Coordinates(1, 0) - Coordinates(1, 1);
  B(0, 5) = 0.f;

  B(1, 0) = 0.f;
  B(1, 1) = Coordinates(0, 2) - Coordinates(0, 1);
  B(1, 2) = 0.f;
  B(1, 3) = Coordinates(0, 0) - Coordinates(0, 2);
  B(1, 4) = 0.f;
  B(1, 5) = Coordinates(0, 1) - Coordinates(0, 0);

  B(2, 0) = Coordinates(0, 2) - Coordinates(0, 1);
  B(2, 1) = Coordinates(1, 1) - Coordinates(1, 2);
  B(2, 2) = Coordinates(0, 0) - Coordinates(0, 2);
  B(2, 3) = Coordinates(1, 2) - Coordinates(1, 0);
  B(2, 4) = Coordinates(0, 1) - Coordinates(0, 0);
  B(2, 5) = Coordinates(1, 0) - Coordinates(1, 1);

  B.divideElementwise(std::abs(area));
}

void Element::calculateCoordinates(CPU_Matrix &Coordinates,
                                   std::vector<Matrix> &nodes) {
  for (size_t i = 0; i < Coordinates.get_numRows(); ++i) {
    for (size_t j = 0; j < Coordinates.get_numCols(); ++j) {
      Coordinates(i, j) = nodes[i][nodesIds[j]];
    }
  }
}

// CalculateStiffnessMatrix will be deprecated!
/*
// 2D only
void Element::CalculateKlocal(Matrix& D, std::vector<Matrix> &nodes) {

  CPU_Matrix Coordinates(DIM, 3);
  calculateCoordinates(Coordinates, nodes);
  float area = calculateArea(Coordinates);
  calculateGradientMatrix(Coordinates, area);


//  Matrix temp1(6, 3);
//  Matrix temp_B(3, 6);


  //K = B.transpose() * D * B * std::abs(C.det()) / 2.0;

//  temp_B = B;

//  temp_B.transpose2();

//  temp1 = temp_B.Product(D);

//  Klocal = temp1.Product(B);

//  Klocal.scale(std::abs(determinant) * 0.5);
}

// 2D only
void Element::CalculateMlocal(float rho, std::vector<Matrix> &nodes, bool lumped) {
  // Zienkiewicz, Taylor, Zhu "The Finite Element Method: Its Basis and Fundamentals" 16.2.3
  float X0 = nodes[0][nodesIds[0]], X1 = nodes[0][nodesIds[1]], X2 = nodes[0][nodesIds[2]];
  float Y0 = nodes[1][nodesIds[0]], Y1 = nodes[1][nodesIds[1]], Y2 = nodes[1][nodesIds[2]];
  float area = 0.5f * std::abs( (X0 - X2) * (Y1 - Y2) - (X1 - X2) * (Y0 - Y2) );
  float mass = rho * area;

  if (lumped) {
    Mlocal(0, 0) = Mlocal(1, 1) = Mlocal(2, 2) = Mlocal(3, 3) = Mlocal(4, 4) = Mlocal(5, 5) = mass / 3;
  } else {
    Mlocal(0, 0) = Mlocal(1, 1) = Mlocal(2, 2) = Mlocal(3, 3) = Mlocal(4, 4) = Mlocal(5, 5) = mass / 6;
    Mlocal(2, 0) = Mlocal(3, 1) = Mlocal(4, 0) = Mlocal(5, 1) = Mlocal(4, 2) = Mlocal(5, 3) =
        Mlocal(0, 2) = Mlocal(1, 3) = Mlocal(0, 4) = Mlocal(1, 5) = Mlocal(2, 4) = Mlocal(3, 5) = mass / 12;
  }
}
*/
void Element::CalculateClocal(float alpha, float beta) {
//  Mlocal.addWeighted(Klocal, alpha, beta);
//  Mlocal.copy(Clocal);
}

// 2D only
void Element::CalculateFlocal2D(BoundaryEdge& edge, std::vector<CPU_Matrix> &nodes, float t) {
  //CheckRunTime(__func__)
  edge.update(t);
  float pressure_value = edge.value;
  float X0 = nodes[0][edge.node[0]], X1 = nodes[0][edge.node[1]];
  float Y0 = nodes[1][edge.node[0]], Y1 = nodes[1][edge.node[1]];
  float edge_length = std::sqrt((X1 - X0) * (X1 - X0) + (Y1 - Y0) * (Y1 - Y0));

  // Grisha: Consider adding the third node when parsing to avoid this abundance of if-else statements
  std::vector<int> edge2elem_num(DIM + 1);
  if (edge.node[0] == nodesIds[0]) {
    edge2elem_num[0] = 0;
  } else if (edge.node[0] == nodesIds[1]) {
    edge2elem_num[0] = 1;
  } else if (edge.node[0] == nodesIds[2]) {
    edge2elem_num[0] = 2;
  }
  if (edge.node[1] == nodesIds[0]) {
    edge2elem_num[1] = 0;
  } else if (edge.node[1] == nodesIds[1]) {
    edge2elem_num[1] = 1;
  } else if (edge.node[1] == nodesIds[2]) {
    edge2elem_num[1] = 2;
  }
  if (edge2elem_num[0] == 0) {
    if (edge2elem_num[1] == 1)
      edge2elem_num[2] = 2;
    else
      edge2elem_num[2] = 1;
  } else if (edge2elem_num[0] == 1) {
    if (edge2elem_num[1] == 0)
      edge2elem_num[2] = 2;
    else
      edge2elem_num[2] = 0;
  } else if (edge2elem_num[0] == 2) {
    if (edge2elem_num[1] == 1)
      edge2elem_num[2] = 0;
    else
      edge2elem_num[2] = 1;
  }

  if (edge.type[0] & Constraint::UX) {
    Flocal[2 * edge2elem_num[0] + 0] = 0.0f;
  } else {
    Flocal[2 * edge2elem_num[0] + 0] = - 0.5f * pressure_value * edge_length * edge.normal[0];
  }
  if (edge.type[0] & Constraint::UY) {
    Flocal[2 * edge2elem_num[0] + 1] = 0.0f;
  } else {
    Flocal[2 * edge2elem_num[0] + 1] = - 0.5f * pressure_value * edge_length * edge.normal[1];
  }

  if (edge.type[1] & Constraint::UX) {
    Flocal[2 * edge2elem_num[1] + 0] = 0.0f;
  } else {
    Flocal[2 * edge2elem_num[1] + 0] = - 0.5f * pressure_value * edge_length * edge.normal[0];
  }
  if (edge.type[1] & Constraint::UY) {
    Flocal[2 * edge2elem_num[1] + 1] = 0.0f;
  } else {
    Flocal[2 * edge2elem_num[1] + 1] = - 0.5f * pressure_value * edge_length * edge.normal[1];
  }

}


void Element::CalculateFlocal3D(BoundaryEdge& edge, std::vector<CPU_Matrix> &nodes, float t) {
  //CheckRunTime(__func__)
  edge.update(t);
  float pressure_value = edge.value;

  float X0 = nodes[0][edge.node[0]], X1 = nodes[0][edge.node[1]], X2 = nodes[0][edge.node[2]];
  float Y0 = nodes[1][edge.node[0]], Y1 = nodes[1][edge.node[1]], Y2 = nodes[1][edge.node[2]];
  float Z0 = nodes[2][edge.node[0]], Z1 = nodes[2][edge.node[1]], Z2 = nodes[2][edge.node[2]];
  struct MathVector{float x, y, z;};
  MathVector a{X1 - X0, Y1 - Y0, Z1 - Z0}, b{X2 - X0, Y2 - Y0, Z2 - Z0};
  float face_area = 0.5f * ((a.y * b.z - a.z * b.y) * (a.y * b.z - a.z * b.y) +
                            (a.z * b.x - a.x * b.z) * (a.z * b.x - a.x * b.z) +
                            (a.x * b.y - a.y * b.x) * (a.x * b.y - a.y * b.x));

  // Grisha: Consider adding the third node when parsing to avoid this abundance of if-else statements
  std::vector<int> edge2elem_num(DIM + 1);
  if (edge.node[0] == nodesIds[0]) {
    edge2elem_num[0] = 0;
  } else if (edge.node[0] == nodesIds[1]) {
    edge2elem_num[0] = 1;
  } else if (edge.node[0] == nodesIds[2]) {
    edge2elem_num[0] = 2;
  } else if (edge.node[0] == nodesIds[3]) {
    edge2elem_num[0] = 3;
  }

  if (edge.node[1] == nodesIds[0]) {
    edge2elem_num[1] = 0;
  } else if (edge.node[1] == nodesIds[1]) {
    edge2elem_num[1] = 1;
  } else if (edge.node[1] == nodesIds[2]) {
    edge2elem_num[1] = 2;
  } else if (edge.node[1] == nodesIds[3]) {
    edge2elem_num[1] = 3;
  }

  if (edge.node[2] == nodesIds[0]) {
    edge2elem_num[2] = 0;
  } else if (edge.node[2] == nodesIds[1]) {
    edge2elem_num[2] = 1;
  } else if (edge.node[2] == nodesIds[2]) {
    edge2elem_num[2] = 2;
  } else if (edge.node[2] == nodesIds[3]) {
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

  if (edge.type[0] & Constraint::UX) {
    Flocal[3 * edge2elem_num[0] + 0] = 0.0f;
  } else {
    Flocal[3 * edge2elem_num[0] + 0] = - 0.5 * pressure_value * face_area * edge.normal[0];
  }
  if (edge.type[0] & Constraint::UY) {
    Flocal[3 * edge2elem_num[0] + 1] = 0.0f;
  } else {
    Flocal[3 * edge2elem_num[0] + 1] = - 0.5 * pressure_value * face_area * edge.normal[1];
  }
  if (edge.type[0] & Constraint::UZ) {
    Flocal[3 * edge2elem_num[0] + 2] = 0.0f;
  } else {
    Flocal[3 * edge2elem_num[0] + 2] = - 0.5 * pressure_value * face_area * edge.normal[2];
  }

  if (edge.type[1] & Constraint::UX) {
    Flocal[3 * edge2elem_num[1] + 0] = 0.0f;
  } else {
    Flocal[3 * edge2elem_num[1] + 0] = - 0.5 * pressure_value * face_area * edge.normal[0];
  }
  if (edge.type[1] & Constraint::UY) {
    Flocal[3 * edge2elem_num[1] + 1] = 0.0f;
  } else {
    Flocal[3 * edge2elem_num[1] + 1] = - 0.5 * pressure_value * face_area * edge.normal[1];
  }
  if (edge.type[1] & Constraint::UZ) {
    Flocal[3 * edge2elem_num[1] + 2] = 0.0f;
  } else {
    Flocal[3 * edge2elem_num[1] + 2] = - 0.5 * pressure_value * face_area * edge.normal[2];
  }

  if (edge.type[2] & Constraint::UX) {
    Flocal[3 * edge2elem_num[2] + 0] = 0.0f;
  } else {
    Flocal[3 * edge2elem_num[2] + 0] = - 0.5 * pressure_value * face_area * edge.normal[0];
  }
  if (edge.type[2] & Constraint::UY) {
    Flocal[3 * edge2elem_num[2] + 1] = 0.0f;
  } else {
    Flocal[3 * edge2elem_num[2] + 1] = - 0.5 * pressure_value * face_area * edge.normal[1];
  }
  if (edge.type[2] & Constraint::UZ) {
    Flocal[3 * edge2elem_num[2] + 2] = 0.0f;
  } else {
    Flocal[3 * edge2elem_num[2] + 2] = - 0.5 * pressure_value * face_area * edge.normal[2];
  }

}


int Element::Global2LocalNode(int glob_num) {
  if (glob_num == this->nodesIds[0])
    return 0;
  if (glob_num == this->nodesIds[1])
    return 1;
  if (glob_num == this->nodesIds[2])
    return 2;
  if (this->DIM == 3 && glob_num == this->nodesIds[3])
    return 3;
  assert(false);
}

void Load::assignElement(int DIM, std::unordered_map <int, std::vector<int>> &nodeAdjElem) {
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

    std::string freq_str = str.substr(0, str.find(","));
    str.erase(0, str.find(",") + 1);

    std::string timeshift_str = str.substr(0, str.find(","));
    str.erase(0, str.find(",") + 1);

    std::string ampl_str = str.substr(0, str.find(")"));

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

    this->timeshift = std::stof(timeshift_str);
    this->freq      = std::stof(freq_str);
  }
}
