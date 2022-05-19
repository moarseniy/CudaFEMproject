
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
  //    B(0, 0) = y[2] - y[3]; B(0, 2) = y[3] - y[1]; B(0, 4) = y[1] - y[2];
  //    B(1, 1) = x[3] - x[2]; B(1, 3) = x[1] - x[3]; B(1, 5) = x[2] - x[1];

  //    B(2, 0) = x[3] - y[2]; B(2, 1) = y[2] - y[3]; B(2, 2) = x[1] - x[3];
  //    B(2, 3) = y[3] - y[1]; B(2, 4) = x[2] - x[1]; B(2, 5) = y[1] - y[2];
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
  double determinant = C.det(3);

  for (int i = 0; i < 3; i++) {
    B(0, 2 * i + 0) = IC(1, i);
    B(0, 2 * i + 1) = 0.0;
    B(1, 2 * i + 0) = 0.0;
    B(1, 2 * i + 1) = IC(2, i);
    B(2, 2 * i + 0) = IC(2, i);
    B(2, 2 * i + 1) = IC(1, i);
  }
  //    B(0, 0) = y[1] - y[2]; B(0, 2) = y[2] - y[0]; B(0, 4) = y[0] - y[1];
  //    B(1, 1) = x[2] - x[1]; B(1, 3) = x[0] - x[2]; B(1, 5) = x[1] - x[0];

  //    B(2, 0) = x[2] - x[1]; B(2, 1) = y[1] - y[2]; B(2, 2) = x[0] - x[2];
  //    B(2, 3) = y[2] - y[0]; B(2, 4) = x[1] - x[0]; B(2, 5) = y[0] - y[1];

  //    float s = std::abs(1.0 / (determinant));
  //    B.scale(s);

  Matrix temp1(6, 3);
  Matrix temp_B(3, 6);


  //K = B.transpose() * D * B * std::abs(C.det()) / 2.0;

  temp_B = B;

  temp_B.transpose2();

  temp1 = temp_B.Product(D);

  Klocal = temp1.Product(B);

  Klocal.scale(std::abs(determinant) * 0.5);
}

void Element::CalculateMlocal(float rho, MyArray& nodesX, MyArray& nodesY, bool lumped) {
  // Zienkiewicz, Taylor, Zhu "The Finite Element Method: Its Basis and Fundamentals" 16.2.3
  float X0 = nodesX[nodesIds[0]], X1 = nodesX[nodesIds[1]], X2 = nodesX[nodesIds[2]];
  float Y0 = nodesY[nodesIds[0]], Y1 = nodesY[nodesIds[1]], Y2 = nodesY[nodesIds[2]];
  float area = 0.5 * std::abs( (X0 - X2) * (Y1 - Y2) - (X1 - X2) * (Y0 - Y2) );
  float mass = rho * area;

  if (lumped) {
    Mlocal(0, 0) = Mlocal(1, 1) = Mlocal(2, 2) = Mlocal(3, 3) = Mlocal(4, 4) = Mlocal(5, 5) = mass / 3;
  } else {
    Mlocal(0, 0) = Mlocal(1, 1) = Mlocal(2, 2) = Mlocal(3, 3) = Mlocal(4, 4) = Mlocal(5, 5) = mass / 6;
    Mlocal(2, 0) = Mlocal(3, 1) = Mlocal(4, 0) = Mlocal(5, 1) = Mlocal(4, 2) = Mlocal(5, 3) =
        Mlocal(0, 2) = Mlocal(1, 3) = Mlocal(0, 4) = Mlocal(1, 5) = Mlocal(2, 4) = Mlocal(3, 5) = mass / 12;
  }
}

void Element::CalculateClocal(float alpha, float beta) {
  if (!(alpha == 0.0f && beta == 0.0f)) Clocal = Mlocal.weightedSum(Klocal, alpha, beta);
}

void Element::CalculateFlocal(BoundaryEdge& edge, MyArray& nodesX, MyArray& nodesY, float t) {
  edge.update(t);
  float pressure_value = edge.value;
  float X0 = nodesX[edge.node0], X1 = nodesX[edge.node1];
  float Y0 = nodesY[edge.node0], Y1 = nodesY[edge.node1];
  float edge_length = std::sqrt((X1 - X0) * (X1 - X0) + (Y1 - Y0) * (Y1 - Y0));

  // Grisha: Consider adding the third node when parsing to avoid this abundance of if-else statements
  int edge2elem_num[3];
  //    std::cout << nodesIds[0] << ' ' << nodesIds[1] << ' ' << nodesIds[2] << ' ' << std::endl;
  //    std::cout << edge.node0 << ' ' << edge.node1 << std::endl;
  if (edge.node0 == nodesIds[0]) {
    edge2elem_num[0] = 0;
  } else if (edge.node0 == nodesIds[1]) {
    edge2elem_num[0] = 1;
  } else if (edge.node0 == nodesIds[2]) {
    edge2elem_num[0] = 2;
  }
  if (edge.node1 == nodesIds[0]) {
    edge2elem_num[1] = 0;
  } else if (edge.node1 == nodesIds[1]) {
    edge2elem_num[1] = 1;
  } else if (edge.node1 == nodesIds[2]) {
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
  } else {  //if (edge2elem_num[0] == 2) {
    if (edge2elem_num[1] == 1)
      edge2elem_num[2] = 0;
    else
      edge2elem_num[2] = 1;
  }

  float X2 = nodesX[nodesIds[edge2elem_num[2]]];
  float Y2 = nodesY[nodesIds[edge2elem_num[2]]];

  //    // Calculate area of element
  //    float area = 0.5 * std::abs( (X0 - X2) * (Y1 - Y2) - (X1 - X2) * (Y0 - Y2) );
  //    float a0 = X1 * Y2 - X2 * Y1, a1 = X2 * Y0 - Y2 * X0;

  //    Flocal[edge2elem_num[0] * 2 + 0] = 0.5 * pressure_value * edge_length * edge.normal_x * (1 +  0.5*a0/area);
  //    Flocal[edge2elem_num[0] * 2 + 1] = 0.5 * pressure_value * edge_length * edge.normal_y * (1 +  0.5*a0/area);
  //    Flocal[edge2elem_num[1] * 2 + 0] = 0.5 * pressure_value * edge_length * edge.normal_x * (1 +  0.5*a1/area);
  //    Flocal[edge2elem_num[1] * 2 + 1] = 0.5 * pressure_value * edge_length * edge.normal_y * (1 +  0.5*a1/area);

  Flocal[2 * edge2elem_num[0] + 0] = - 0.5 * pressure_value * edge_length * edge.normal_x;
  Flocal[2 * edge2elem_num[0] + 1] = - 0.5 * pressure_value * edge_length * edge.normal_y;
  Flocal[2 * edge2elem_num[1] + 0] = - 0.5 * pressure_value * edge_length * edge.normal_x;
  Flocal[2 * edge2elem_num[1] + 1] = - 0.5 * pressure_value * edge_length * edge.normal_y;

}

int Element::Global2LocalNode(int glob_num) {
  if (glob_num == this->nodesIds[0])
    return 0;
  if (glob_num == this->nodesIds[1])
    return 1;
  if (glob_num == this->nodesIds[2])
    return 2;
  assert(false);
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
        std::cos(w0*t) * std::sqrtf(3.0f) * (t*t/w1 + t/(w1*w1)) );

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
    std::cout << "DEBUG : xstr = " << str << std::endl;
    std::cout << "DEBUG : wavelet_str = " << wavelet_str << std::endl;

    std::string freq_str = str.substr(0, str.find(","));
    str.erase(0, str.find(",") + 1);
    std::cout << "DEBUG : xstr = " << str << std::endl;
    std::cout << "DEBUG : freq_str = " << freq_str << std::endl;

    std::string timeshift_str = str.substr(0, str.find(","));
    str.erase(0, str.find(",") + 1);
    std::cout << "DEBUG : xstr = " << str << std::endl;
    std::cout << "DEBUG : timeshift_str = " << timeshift_str << std::endl;

    std::string ampl_str = str.substr(0, str.find(")"));
    //xstr.erase(0, xstr.find(")") + 1);
    std::cout << "DEBUG : xstr = " << str << std::endl;
    std::cout << "DEBUG : ampl_str = " << ampl_str << std::endl;

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

    std::cout << "DEBUG : wavelet_str = " << wavelet_str << std::endl;

    this->timeshift = std::stof(timeshift_str);
    this->freq      = std::stof(freq_str);
    std::cout << "value = " << this->value << "; timeshift = " << this->timeshift << "; load.freq = " << this->freq << std::endl;

  }
}
