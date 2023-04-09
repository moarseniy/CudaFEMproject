#include <datakeeper.h>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

void dataKeeper::parseJsonConfig(std::string configPath) {
  std::ifstream ifs(configPath);
  json config = json::parse(ifs);

  DIM = config.at("DIM");

  paths.proj_dir = config.value("CudaFEMproject_path", "undefined_path");
  paths.work_dir = config.value("working_path", "undefined_path");

  paths.res_dir = fs::joinString(paths.work_dir, "/fem_results/" + std::to_string(DIM) + "D/" + taskName);
  paths.prep_mesh_dir = fs::joinString(paths.work_dir, "/prepared_meshes/" + std::to_string(DIM) + "D/" + taskName);
  paths.output_vtk = fs::joinString(paths.res_dir, "/results.vtk");

  // required for static problems
  mechParams.youngModulus = config.value("youngModulus", 2e+7f);
  mechParams.poissonRatio = config.value("poissonRatio", 0.3f);

  // required for dynamic problems
  taskType = config.value("task_type", "static");
  if (taskType == "dynamic") {
    mechParams.rho = config.value("rho", 2400.f);
    mechParams.beta1 = config.value("beta1", 0.5f);
    mechParams.beta2 = config.value("beta2", 0.f);
    mechParams.damping_alpha = config.value("damping_alpha", 0.f);
    mechParams.damping_beta = config.value("damping_beta", 0.f);
    mechParams.dt = config.value("dt", 1e-3f);
    mechParams.endtime = config.value("endtime", 0.1f);
  }
}

void dataKeeper::parseCounters() {
  nodes_file.open         (fs::joinString(paths.prep_mesh_dir, "/nodes.txt"),         std::fstream::in);
  elements_file.open      (fs::joinString(paths.prep_mesh_dir, "/elements.txt"),      std::fstream::in);
  loads_file.open         (fs::joinString(paths.prep_mesh_dir, "/loads.txt"),         std::fstream::in);
  constraints_file.open   (fs::joinString(paths.prep_mesh_dir, "/constraints.txt"),   std::fstream::in);
  stress_file.open        (fs::joinString(paths.prep_mesh_dir, "/stress.txt"),        std::fstream::in);

  nodes_file          >> nodesCount;
  elements_file       >> elementsCount;
  constraints_file    >> constraintsCount;
  loads_file          >> loadsCount;
  stress_file         >> boundaryEdgesCount;
}

void dataKeeper::allocateMemory() {
//  nodes->resize(nodesCount, DIM);
//  elementsIds.resize(elementsCount, DIM + 1);
//  D.resize(3 * (DIM - 1), 3 * (DIM - 1));
//  constraintsIds.resize(1, constraintsCount);
  nodes = Matrix::setMatrix(CPU, nodesCount, DIM);
  elementsIds = Matrix::setMatrix(CPU, elementsCount, DIM + 1);
  D = Matrix::setMatrix(CPU, 3 * (DIM - 1), 3 * (DIM - 1));

//  constraintsIds = Matrix::setMatrix(CPU, nodesCount, DIM);
//  constraintsIds->setTo(0);
//  constraintsIds = Matrix::setMatrix(CPU);
  constraintsTypes = Matrix::setVector(CPU, nodesCount);
  constraintsTypes->setTo(0);

  boundaryNormals = Matrix::setMatrix(CPU, boundaryEdgesCount, DIM);
  boundaryAdjElems = Matrix::setMatrix(CPU, 1, boundaryEdgesCount);
  boundaryNodes = Matrix::setMatrix(CPU, boundaryEdgesCount, DIM + 1);
  boundaryPressureValues = Matrix::setMatrix(CPU, 1, boundaryEdgesCount);

  displacements = Matrix::setMatrix(CPU, DIM, nodesCount);
//  pressure = Matrix::setVector(CPU, boundaryEdgesCount);
}

dataKeeper::~dataKeeper() {
  if (nodes)
    delete nodes;
  if (elementsIds)
    delete elementsIds;
  if (D)
    delete D;
  if (constraintsIds)
    delete constraintsIds;
  if (constraintsTypes)
    delete constraintsTypes;

  if (boundaryNormals)
    delete boundaryNormals;
  if (boundaryAdjElems)
    delete boundaryAdjElems;
  if (boundaryNodes)
    delete boundaryNodes;
  if (boundaryPressureValues)
    delete boundaryPressureValues;
  if (displacements)
    delete displacements;
}
void dataKeeper::CreateMatrixD() {
  float poissonRatio = mechParams.poissonRatio;
  float youngModulus = mechParams.youngModulus;
  if (DIM == 2) {
    (*D)(0, 0) = 1.0;          (*D)(0, 1) = poissonRatio; (*D)(0, 2) = 0.0;
    (*D)(1, 0) = poissonRatio; (*D)(1, 1) = 1.0;          (*D)(1, 2) = 0.0;
    (*D)(2, 0) = 0.0;          (*D)(2, 1) = 0.0;          (*D)(2, 2) = (1.0 - poissonRatio) / 2.0;
    D->scale(youngModulus / (1.0 - pow(poissonRatio, 2.0)));
  } else if (DIM == 3) {
    (*D)(0, 0) = (*D)(1, 1) = (*D)(2, 2) = 1.0;
    (*D)(0, 1) = (*D)(1, 0) = (*D)(0, 2) = (*D)(2, 0) = (*D)(2, 1) = (*D)(1, 2) = poissonRatio / (1.0 - poissonRatio);
    (*D)(3, 3) = (*D)(4, 4) = (*D)(5, 5) = (1.0 - 2.0 * poissonRatio) / (2.0 * (1.0 - poissonRatio));
    D->scale((youngModulus * (1.0 - poissonRatio)) / ((1.0 + poissonRatio) * (1.0 - 2.0 * poissonRatio)));
  } else {
    throw std::runtime_error("Error in FEMdataKeeper::CreateMatrixD()");
  }
}

void dataKeeper::ParseNodes() {
  CheckRunTime(__func__)
  for (size_t i = 0; i < nodes->get_numElements(); ++i) {
      nodes_file >> (*nodes)[i];
  }
}

void dataKeeper::ParseElements() {
  CheckRunTime(__func__)
  for (size_t i = 0; i < elementsIds->get_numElements(); ++i) {
    this->elements_file >> (*elementsIds)[i];
  }

// TODO: Jacobian calculation
//        float X2 = nodesX[element.nodesIds[1]] - nodesX[element.nodesIds[0]];
//        float X3 = nodesX[element.nodesIds[2]] - nodesX[element.nodesIds[0]];
//        float Y2 = nodesY[element.nodesIds[1]] - nodesY[element.nodesIds[0]];
//        float Y3 = nodesY[element.nodesIds[2]] - nodesY[element.nodesIds[0]];
//        element.jacobian = 1.0 / (X3 * Y2 - X2 * Y3);
}

void dataKeeper::ParseConstraints() {
  CheckRunTime(__func__)

  int node, type;
  Constraint::Type constraintType;
  std::vector<float> temp;
  for (size_t i = 0; i < constraintsCount; ++i) {
    constraints_file >> node >> type;
    (*constraintsTypes)[node] = type;

    constraintType = static_cast<Constraint::Type>(type);
    if (constraintType & Constraint::UX) {
      temp.push_back(DIM * node + 0);
    }
    if (constraintType & Constraint::UY) {
      temp.push_back(DIM * node + 1);
    }
    if (constraintType & Constraint::UZ) {
      temp.push_back(DIM * node + 2);
    }
  }

  if (temp.size() > 0) {
    constraintsIds = new CPU_Matrix(temp.size(), 1);
    std::copy(temp.data(), temp.data() + temp.size(), constraintsIds->get_data());
  }
  std::cout << "Constraints ID's (" << temp.size() << ") found\n";
}

void dataKeeper::ParseLoads() {
  CheckRunTime(__func__)
  for (int i = 0; i < this->loadsCount; ++i) {
    int node;
    std::string xstr, ystr, str;
    std::string zstr; // if 3D

    this->loads_file >> node;
    std::getline(this->loads_file, str);
    // https://stackoverflow.com/a/39080627
    std::regex wavelet_expression("(\\w+(\\((.*?)\\))|[+-]?(\\d+([.]\\d*)?([eE][+-]?\\d+)?|[.]\\d+([eE][+-]?\\d+)?))");
    std::smatch res;

    int loop_it = 0;
    while (std::regex_search(str, res, wavelet_expression)) {
      ++loop_it;
      if (loop_it == 1)      xstr = res[0];
      else if (loop_it == 2) ystr = res[0];
      else if (this->DIM == 3 && loop_it == 3) zstr = res[0];
      str = res.suffix();
    }
    // X
    Load loadX;
    loadX.parseString(xstr);
    if (!(loadX.ampl == 0.0f)) {
      loadX.dof = node * this->DIM + 0;
      this->loads.push_back(loadX);
    }
    // Y
    Load loadY;
    loadY.parseString(ystr);
    if (!(loadY.ampl == 0.0f)) {
      loadY.dof = DIM == node * this->DIM + 1;
      this->loads.push_back(loadY);
    }

    // Z
    if (this->DIM == 3) {
      Load loadZ;
      loadZ.parseString(zstr);
      if (!(loadZ.ampl == 0.0f)) {
        loadZ.dof = node * DIM + 2;
        this->loads.push_back(loadZ);
      }
    }

  }
}

void dataKeeper::ParseBoundaryEdges() {
  CheckRunTime(__func__)
  for (size_t i = 0; i < boundaryEdgesCount; ++i) {
    for (size_t j = 0; j < DIM + 1; ++j) {
      stress_file >> (*boundaryNodes)(i, j);
    }
    stress_file >> (*boundaryAdjElems)[i];

    std::vector<float> normal_vec(DIM);
    float normal_length = 0.f;
    for (size_t j = 0; j < DIM; ++j) {
      stress_file >> normal_vec[j];
//      std::cout << normal_vec[j] << " ";
      normal_length += normal_vec[j] * normal_vec[j];
    }
    for (size_t j = 0; j < DIM; ++j) {
      (*boundaryNormals)(i, j) = normal_vec[j] / std::sqrt(normal_length);
    }

    stress_file >> (*boundaryPressureValues)[i];
  }
//  boundaryNormals->Show();


  /*
  for (size_t i = 0; i < boundaryEdgesCount; ++i) {
    BoundaryEdge edge(DIM);
    std::vector<float> normal_vec(DIM);

    for (size_t j = 0; j < DIM; ++j) {
      stress_file >> edge.node[j];
    }
    stress_file >> edge.adj_elem1;
    for (size_t j = 0; j < DIM; ++j) {
      stress_file >> normal_vec[j];
    }
    for (size_t j = 0; j < constraintsCount; ++j) {
      for (size_t k = 0; k < DIM; ++k) {
        if (constraints[j].node == edge.node[k]) {
          edge.type[k] = constraints[j].type;
        }
      }
    }

    std::string str;
    std::getline(stress_file, str);
    // https://stackoverflow.com/a/39080627
    std::regex wavelet_expression("(\\w+(\\((.*?)\\))|[+-]?(\\d+([.]\\d*)?([eE][+-]?\\d+)?|[.]\\d+([eE][+-]?\\d+)?))");
    std::smatch res;
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

    if (!(edge.ampl == 0.0f))
      this->boundary.push_back(edge);
  }
  */
}

void dataKeeper::parseFiles() {
CheckRunTime(__func__)
  CreateMatrixD();
  ParseNodes();
  ParseElements();
  ParseConstraints();
  ParseLoads();
  ParseBoundaryEdges();
}

void dataKeeper::ShowInfo() {
  std::cout << "==========INFO==========" <<
               "\nTask name = " << taskName <<
               "\nNodes count = " << nodesCount <<
               "\nElements count = " << elementsCount <<
               "\nConstraints count = " << constraintsCount <<
               "\nLoads count = " << loadsCount <<
               "\nBoundary edges with b.c. count = " << boundaryEdgesCount <<
               "\n========================\n\n";
}

dataKeeper::dataKeeper(std::string configPath, std::string task) {
  taskName = task;
  parseJsonConfig(configPath);
  parseCounters();
  allocateMemory();
  parseFiles();
}

FEMdataKeeper::FEMdataKeeper(std::string name,
                             int DIM,
                             std::string project_directory,
                             std::string prepared_meshes_directory,
                             std::string results_directory,
                             float poissonRatio, float youngModulus) {
  this->name = name;
  this->DIM = DIM;
  this->proj_dir = project_directory;
  this->prep_mesh_dir = prepared_meshes_directory;
  this->res_dir = results_directory;
  this->poissonRatio = poissonRatio; this->youngModulus = youngModulus;

  this->nodes_file.open         (this->prep_mesh_dir + "/nodes.txt",         std::fstream::in);
  this->elements_file.open      (this->prep_mesh_dir + "/elements.txt",      std::fstream::in);
  this->loads_file.open         (this->prep_mesh_dir + "/loads.txt",         std::fstream::in);
  this->constraints_file.open   (this->prep_mesh_dir + "/constraints.txt",   std::fstream::in);
  this->stress_file.open        (this->prep_mesh_dir + "/stress.txt",        std::fstream::in);

  this->nodes_file          >> this->nodesCount;
  this->elements_file       >> this->elementsCount;
  this->constraints_file    >> this->constraintsCount;
  this->loads_file          >> this->loadsCount;
  this->stress_file         >> this->boundaryEdgesCount;

  // Allocate dynamic memory
  for (int i = 0; i < this->DIM; ++i) {
    this->nodes.push_back(CPU_Matrix(this->nodesCount, 1));
  }

  this->D.resize(3 * (DIM - 1), 3 * (DIM - 1));
  this->displacements = Matrix::setMatrix(CPU, this->DIM, this->nodesCount);
  this->pressure->setVector(CPU, this->boundaryEdgesCount);
  this->all_B.resize(3 * (this->DIM - 1) * 6 * (this->DIM - 1) * this->elementsCount, 1);

  // Parse the files
  this->parseFiles();
}

void FEMdataKeeper::ShowInfo() {
  std::cout << "==========INFO==========" <<
               "\nTask name = " << name <<
               "\nNodes count = " << nodesCount <<
               "\nElements count = " << elementsCount <<
               "\nConstraints count = " << constraintsCount <<
               "\nLoads count = " << loadsCount <<
               "\nBoundary edges with b.c. count = " << boundaryEdgesCount <<
               "\n========================\n\n";
}

void FEMdataKeeper::ParseNodes() {
  CheckRunTime(__func__)
  for (int i = 0; i < this->nodesCount; ++i) {
    for (int j = 0; j < this->DIM; ++j) {
      this->nodes_file >> this->nodes[j][i];
    }
  }
}

void FEMdataKeeper::ParseElements() {
  CheckRunTime(__func__)
  for (int i = 0; i < this->elementsCount; ++i) {
    Element element(this->DIM);
    for (int j = 0; j < this->DIM + 1; ++j) {
      int nodeId;
      this->elements_file >> nodeId;
      element.nodesIds[j] = nodeId;
    }

    // ToDO in future: Jacobian calculation
    //        float X2 = nodesX[element.nodesIds[1]] - nodesX[element.nodesIds[0]];
    //        float X3 = nodesX[element.nodesIds[2]] - nodesX[element.nodesIds[0]];
    //        float Y2 = nodesY[element.nodesIds[1]] - nodesY[element.nodesIds[0]];
    //        float Y3 = nodesY[element.nodesIds[2]] - nodesY[element.nodesIds[0]];
    //        element.jacobian = 1.0 / (X3 * Y2 - X2 * Y3);
    this->elements.push_back(element);
  }
}

void FEMdataKeeper::ParseConstraints() {
  CheckRunTime(__func__)
  for (int i = 0; i < this->constraintsCount; ++i) {
    Constraint constraint;
    int type;
    this->constraints_file >> constraint.node >> type;
    constraint.type = static_cast<Constraint::Type>(type);
    this->constraints.push_back(constraint);
  }
}

void FEMdataKeeper::ParseLoads() {
  CheckRunTime(__func__)
  for (int i = 0; i < this->loadsCount; ++i) {
    int node;
    std::string xstr, ystr, str;
    std::string zstr; // if 3D

    this->loads_file >> node;
    std::getline(this->loads_file, str);
    // https://stackoverflow.com/a/39080627
    std::regex wavelet_expression("(\\w+(\\((.*?)\\))|[+-]?(\\d+([.]\\d*)?([eE][+-]?\\d+)?|[.]\\d+([eE][+-]?\\d+)?))");
    std::smatch res;

    int loop_it = 0;
    while (std::regex_search(str, res, wavelet_expression)) {
      ++loop_it;
      if (loop_it == 1)      xstr = res[0];
      else if (loop_it == 2) ystr = res[0];
      else if (this->DIM == 3 && loop_it == 3) zstr = res[0];
      str = res.suffix();
    }
    // X
    Load loadX;
    loadX.parseString(xstr);
    if (!(loadX.ampl == 0.0f)) {
      loadX.dof = node * this->DIM + 0;
      this->loads.push_back(loadX);
    }
    // Y
    Load loadY;
    loadY.parseString(ystr);
    if (!(loadY.ampl == 0.0f)) {
      loadY.dof = DIM == node * this->DIM + 1;
      this->loads.push_back(loadY);
    }

    // Z
    if (this->DIM == 3) {
      Load loadZ;
      loadZ.parseString(zstr);
      if (!(loadZ.ampl == 0.0f)) {
        loadZ.dof = node * DIM + 2;
        this->loads.push_back(loadZ);
      }
    }

  }
}

void FEMdataKeeper::ParseBoundaryEdges() {
  CheckRunTime(__func__)
  for (int i = 0; i < this->boundaryEdgesCount; ++i) {
    BoundaryEdge edge(this->DIM);
    std::vector<float> normal_vec(this->DIM);

    for (int j = 0; j < this->DIM; ++j) {
      this->stress_file >> edge.node[j];
    }
    int temp;
    stress_file >> temp;

    this->stress_file >> edge.adj_elem1;
    for (int j = 0; j < this->DIM; ++j) {
      this->stress_file >> normal_vec[j];
    }

    for (std::vector<Constraint>::const_iterator it = this->constraints.begin();
         it != this->constraints.end();
         ++it) {
      for (int j = 0; j < this->DIM; ++j) {
        if (it->node == edge.node[j]) {
          edge.type[j] = it->type;
        }
      }
    }

    std::string str;
    std::getline(stress_file, str);
    // https://stackoverflow.com/a/39080627
    std::regex wavelet_expression("(\\w+(\\((.*?)\\))|[+-]?(\\d+([.]\\d*)?([eE][+-]?\\d+)?|[.]\\d+([eE][+-]?\\d+)?))");
    std::smatch res;
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

    if (!(edge.ampl == 0.0f)) this->boundary.push_back(edge);
  }
}

void FEMdataKeeper::CreateMatrixD() {
  if (DIM == 2) {
    D(0, 0) = 1.0;          D(0, 1) = poissonRatio; D(0, 2) = 0.0;
    D(1, 0) = poissonRatio;	D(1, 1) = 1.0;          D(1, 2) = 0.0;
    D(2, 0) = 0.0;          D(2, 1) = 0.0;          D(2, 2) = (1.0 - poissonRatio) / 2.0;
    D.scale(youngModulus / (1.0 - pow(poissonRatio, 2.0)));
  } else if (DIM == 3) {
    D(0, 0) = D(1, 1) = D(2, 2) = 1.0;
    D(0, 1) = D(1, 0) = D(0, 2) = D(2, 0) = D(2, 1) = D(1, 2) = poissonRatio / (1.0 - poissonRatio);
    D(3, 3) = D(4, 4) = D(5, 5) = (1.0 - 2.0 * poissonRatio) / (2.0 * (1.0 - poissonRatio));
    D.scale((youngModulus * (1.0 - poissonRatio)) / ((1.0 + poissonRatio) * (1.0 - 2.0 * poissonRatio)));
  } else {
    throw std::runtime_error("Error in FEMdataKeeper::CreateMatrixD()");
  }
}

void FEMdataKeeper::parseFiles() {
  CheckRunTime(__func__)
  this->CreateMatrixD();
  this->ParseNodes();
  this->ParseElements();
  this->ParseConstraints();
  this->ParseLoads();
  this->ParseBoundaryEdges();
}


//void ResultsDataKeeper::MakeResults() {
//  this->CalculateStressAndDeformation();

//  float fixed_value = 1.0;// x -3.0; y -3.0
//  float a = 0.0f;// x -3.0 y -4.0
//  float b = 8.0f;// x 5.0 y -3.0;

//  if (this->withStressAlongAxis) {
//    CalculateStressAlongAxis(this->StressComponents,
//                           "x", "xy", fixed_value, a, b,
//                           this->Stress, this->FEMdata->nodes, this->FEMdata->elements);
//  }
//  if (this->withSmooth) {
//    SmoothResults("xx", this->SmoothStress, this->Stress,
//                  this->FEMdata->nodesCount, this->FEMdata->nodes, this->FEMdata->elements);
//    CalculateStressAlongAxisSmooth(this->StressComponentsSmooth,
//                                   "x", fixed_value, a, b, this->SmoothStress,
//                                   this->FEMdata->nodes, this->FEMdata->elements);
//  }

//  a = 0.0f;
//  b = 6.0f;
//  float k = -0.6655f;
//  float m = -0.0035f;

//  if (this->withMises) {
//    CalculateMisesAlongLineMises(this->MisesComponents,
//                                 k, m, a, b, this->sigma_mises,
//                                 this->FEMdata->nodes, this->FEMdata->elements);
//  }
//}

//void ResultsDataKeeper::CalculateStressAndDeformation() {
//  CheckRunTime(__func__)

//  int DIM = this->FEMdata->DIM;
//  MyArray StressVector(3 * (DIM - 1));
//  MyArray DeformationVector(3 * (DIM - 1));
//  MyArray delta(6 * (DIM - 1));

//  Matrix B(3 * (DIM - 1), 6 * (DIM - 1));

//  std::vector<Element>  &elements       = this->FEMdata->elements;
//  MyArray               &displacements  = this->FEMdata->displacements;
//  int i = 0;
//  for (auto elem = elements.begin(); elem != elements.end(); ++elem) {
//  //for (int i = 0; i < elements.size(); ++i) {
//    std::vector<int> &nodesIds = elem->nodesIds;
//    if (DIM == 2) {
//      delta[0] = displacements[2 * nodesIds[0] + 0];
//      delta[1] = displacements[2 * nodesIds[0] + 1];
//      delta[2] = displacements[2 * nodesIds[1] + 0];
//      delta[3] = displacements[2 * nodesIds[1] + 1];
//      delta[4] = displacements[2 * nodesIds[2] + 0];
//      delta[5] = displacements[2 * nodesIds[2] + 1];
//    } else if (DIM == 3) {
//      delta[0] = displacements[3 * nodesIds[0] + 0];
//      delta[1] = displacements[3 * nodesIds[0] + 1];
//      delta[2] = displacements[3 * nodesIds[0] + 2];
//      delta[3] = displacements[3 * nodesIds[1] + 0];
//      delta[4] = displacements[3 * nodesIds[1] + 1];
//      delta[5] = displacements[3 * nodesIds[1] + 2];
//      delta[6] = displacements[3 * nodesIds[2] + 0];
//      delta[7] = displacements[3 * nodesIds[2] + 1];
//      delta[8] = displacements[3 * nodesIds[2] + 2];
//      delta[9] = displacements[3 * nodesIds[3] + 0];
//      delta[10] = displacements[3 * nodesIds[3] + 1];
//      delta[11] = displacements[3 * nodesIds[3] + 2];
//    }

//    for (int k = 0; k < 3 * (DIM - 1); ++k)
//      for (int j = 0; j < 6 * (DIM - 1); ++j)
//        // variable i appers only here
//        B(k, j) = this->FEMdata->all_B[j + 6 * (DIM - 1) * k + 3 * (DIM - 1) * 6 * (DIM - 1) * i];

//    //DeformationVector = it->B.Product(delta);
//    DeformationVector = B.Product(delta);
//    StressVector = this->FEMdata->D.Product(DeformationVector);

//    Deformation.push_back(DeformationVector);
//    Stress.push_back(StressVector);

//    ++i;
//  }
//}
