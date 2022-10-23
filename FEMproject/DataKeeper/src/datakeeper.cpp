#include "datakeeper.h"

void FEMdataKeeper::ParseFiles(float poissonRatio, float youngModulus) {
  CheckRunTime(__func__)

  if (DIM == 2) {
    D(0,0) = 1.0;			D(0, 1) = poissonRatio;	D(0, 2) = 0.0;
    D(1, 0) = poissonRatio;	D(1, 1) = 1.0; 			D(1, 2) = 0.0;
    D(2, 0) = 0.0;        	D(2, 1) = 0.0;        	D(2, 2) = (1.0 - poissonRatio) / 2.0;
    D.scale(youngModulus / (1.0 - pow(poissonRatio, 2.0)));
  } else if (DIM == 3) {
    D(0, 0) = D(1, 1) = D(2, 2) = 1.0;
    D(0, 1) = D(1, 0) = D(0, 2) = D(2, 0) = D(2, 1) = D(1, 2) = poissonRatio / (1.0 - poissonRatio);
    D(3, 3) = D(4, 4) = D(5, 5) = (1.0 - 2.0 * poissonRatio) / (2.0 * (1.0 - poissonRatio));
    D.scale((youngModulus * (1.0 - poissonRatio)) / ((1.0 + poissonRatio) * (1.0 - 2.0 * poissonRatio)));
  } else {
    throw std::runtime_error("DIM error");
  }

  // NODES
  for (int i = 0; i < nodesCount; ++i) {
    for (int j = 0; j < DIM; ++j) {
      nodes_file >> nodes[j][i];
    }
  }

  // ELEMENTS
  for (int i = 0; i < elementsCount; ++i) {
    Element element(DIM);
    for (int j = 0; j < DIM + 1; ++j) {
      elements_file >> element.nodesIds[j];
    }

    // ToDO in future: Jacobian calculation
    //        float X2 = nodesX[element.nodesIds[1]] - nodesX[element.nodesIds[0]];
    //        float X3 = nodesX[element.nodesIds[2]] - nodesX[element.nodesIds[0]];
    //        float Y2 = nodesY[element.nodesIds[1]] - nodesY[element.nodesIds[0]];
    //        float Y3 = nodesY[element.nodesIds[2]] - nodesY[element.nodesIds[0]];
    //        element.jacobian = 1.0 / (X3 * Y2 - X2 * Y3);
    elements.push_back(element);
  }

  // CONSTRAINTS
  for (int i = 0; i < constraintsCount; ++i) {
    Constraint constraint;
    int type;
    constraints_file >> constraint.node >> type;
    constraint.type = static_cast<Constraint::Type>(type);
    constraints.push_back(constraint);
  }

  // LOADS
  for (int i = 0; i < loadsCount; ++i) {
    int node;
    std::string xstr, ystr, str;
    std::string zstr; // if 3D

    loads_file >> node;
    std::getline(loads_file, str);
    // https://stackoverflow.com/a/39080627
    std::regex wavelet_expression("(\\w+(\\((.*?)\\))|[+-]?(\\d+([.]\\d*)?([eE][+-]?\\d+)?|[.]\\d+([eE][+-]?\\d+)?))");
    smatch res;

    int loop_it = 0;
    while (std::regex_search(str, res, wavelet_expression)) {
      ++loop_it;
      if (loop_it == 1)      xstr = res[0];
      else if (loop_it == 2) ystr = res[0];
      else if (DIM == 3 && loop_it == 3) zstr = res[0];
      str = res.suffix();
    }
    // X
    Load loadX;
    loadX.parseString(xstr);
    if (!(loadX.ampl == 0.0f)) {
      loadX.dof = node * DIM + 0;
      loads.push_back(loadX);
    }
    // Y
    Load loadY;
    loadY.parseString(ystr);
    if (!(loadY.ampl == 0.0f)) {
      loadY.dof = DIM == node * DIM + 1;
      loads.push_back(loadY);
    }

    // Z
    if (DIM == 3) {
      Load loadZ;
      loadZ.parseString(zstr);
      if (!(loadZ.ampl == 0.0f)) {
        loadZ.dof = node * DIM + 2;
        loads.push_back(loadZ);
      }
    }

  }

  // BOUNDARY EDGES
  for (int i = 0; i < boundaryEdgesCount; ++i) {
    BoundaryEdge edge(DIM);
    MyArray normal_vec(DIM);

    for (int j = 0; j < DIM; ++j) {
      stress_file >> edge.node[j];
    }
    stress_file >> edge.adj_elem1;
    for (int j = 0; j < DIM; ++j) {
      stress_file >> normal_vec[j];
    }

    for (std::vector<Constraint>::const_iterator it = constraints.begin(); it != constraints.end(); ++it) {
      for (int j = 0; j < DIM; ++j) {
        if (it->node == edge.node[j]) {
          edge.type[j] = it->type;
        }
      }
    }

    string str;
    std::getline(stress_file, str);
    // https://stackoverflow.com/a/39080627
    std::regex wavelet_expression("(\\w+(\\((.*?)\\))|[+-]?(\\d+([.]\\d*)?([eE][+-]?\\d+)?|[.]\\d+([eE][+-]?\\d+)?))");
    smatch res;
    std::regex_search(str, res, wavelet_expression);
    str = res[0];
    edge.parseString(str);

    float normal_length = 0.f;
    for (int j = 0; j < DIM; ++j) {
      normal_length += normal_vec[j] * normal_vec[j];
    }
    for (int j = 0; j < DIM; ++j) {
      edge.normal[j] = normal_vec[j] / std::sqrt(normal_length);
    }

    if (!(edge.ampl == 0.0f)) boundary.push_back(edge);
  }
}
