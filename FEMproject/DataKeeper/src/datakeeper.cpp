#include "datakeeper.h"

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

  this->nodes_file.open         (this->prep_mesh_dir + "/nodes.txt",         fstream::in);
  this->elements_file.open      (this->prep_mesh_dir + "/elements.txt",      fstream::in);
  this->loads_file.open         (this->prep_mesh_dir + "/loads.txt",         fstream::in);
  this->constraints_file.open   (this->prep_mesh_dir + "/constraints.txt",   fstream::in);
  this->stress_file.open        (this->prep_mesh_dir + "/stress.txt",        fstream::in);

  this->nodes_file          >> this->nodesCount;
  this->elements_file       >> this->elementsCount;
  this->constraints_file    >> this->constraintsCount;
  this->loads_file          >> this->loadsCount;
  this->stress_file         >> this->boundaryEdgesCount;

  // Allocate dynamic memory
  for (int i = 0; i < this->DIM; ++i) {
    MyArray nodes_vec(this->nodesCount);
    this->nodes.push_back(nodes_vec);
  }
  this->displacements.Resize(this->DIM * this->nodesCount);
  this->D.Resize(3 * (this->DIM - 1), 3 * (this->DIM - 1));
  this->pressure.Resize(this->boundaryEdgesCount);
  this->all_B.Resize(3 * (this->DIM - 1) * 6 * (this->DIM - 1) * this->elementsCount);

  // Parse the files
  this->Parse();
}

void FEMdataKeeper::ShowInfo() {
  std::cout << "==========INFO==========" <<
               "\nTask name = " << this->name <<
               "\nNodes count = " << this->nodesCount <<
               "\nElements count = " << this->elementsCount <<
               "\nConstraints count = " << this->constraintsCount <<
               "\nLoads count = " << this->loadsCount <<
               "\nBoundary edges with b.c. count = " << this->boundaryEdgesCount <<
               "\n========================\n\n";
}

void FEMdataKeeper::ParseNodes()
{
  CheckRunTime(__func__)
  for (int i = 0; i < this->nodesCount; ++i) {
    for (int j = 0; j < this->DIM; ++j) {
      this->nodes_file >> this->nodes[j][i];
    }
  }
}

void FEMdataKeeper::ParseElements()
{
  CheckRunTime(__func__)
  for (int i = 0; i < this->elementsCount; ++i) {
    Element element(this->DIM);
    for (int j = 0; j < this->DIM + 1; ++j) {
      this->elements_file >> element.nodesIds[j];
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

void FEMdataKeeper::ParseConstraints()
{
  CheckRunTime(__func__)
  for (int i = 0; i < this->constraintsCount; ++i) {
    Constraint constraint;
    int type;
    this->constraints_file >> constraint.node >> type;
    constraint.type = static_cast<Constraint::Type>(type);
    this->constraints.push_back(constraint);
  }
}

void FEMdataKeeper::ParseLoads()
{
  CheckRunTime(__func__)
  for (int i = 0; i < this->loadsCount; ++i) {
    int node;
    std::string xstr, ystr, str;
    std::string zstr; // if 3D

    this->loads_file >> node;
    std::getline(this->loads_file, str);
    // https://stackoverflow.com/a/39080627
    std::regex wavelet_expression("(\\w+(\\((.*?)\\))|[+-]?(\\d+([.]\\d*)?([eE][+-]?\\d+)?|[.]\\d+([eE][+-]?\\d+)?))");
    smatch res;

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

void FEMdataKeeper::ParseBoundaryEdges()
{
  CheckRunTime(__func__)
  for (int i = 0; i < this->boundaryEdgesCount; ++i) {
    BoundaryEdge edge(this->DIM);
    MyArray normal_vec(this->DIM);

    for (int j = 0; j < this->DIM; ++j) {
      this->stress_file >> edge.node[j];
    }
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

    if (!(edge.ampl == 0.0f)) this->boundary.push_back(edge);
  }
}

void FEMdataKeeper::CreateMatrixD()
{
  CheckRunTime(__func__)
  if (this->DIM == 2) {
    this->D(0,0) = 1.0;                 this->D(0, 1) = this->poissonRatio;	this->D(0, 2) = 0.0;
    this->D(1, 0) = this->poissonRatio;	this->D(1, 1) = 1.0;                this->D(1, 2) = 0.0;
    this->D(2, 0) = 0.0;                this->D(2, 1) = 0.0;                this->D(2, 2) = (1.0 - this->poissonRatio) / 2.0;
    this->D.scale(this->youngModulus / (1.0 - pow(this->poissonRatio, 2.0)));
  } else if (this->DIM == 3) {
    this->D(0, 0) = this->D(1, 1) = this->D(2, 2) = 1.0;
    this->D(0, 1) = this->D(1, 0) = this->D(0, 2) = this->D(2, 0)
        = this->D(2, 1) = this->D(1, 2) = this->poissonRatio / (1.0 - this->poissonRatio);
    this->D(3, 3) = this->D(4, 4) = this->D(5, 5) = (1.0 - 2.0 * this->poissonRatio) / (2.0 * (1.0 - this->poissonRatio));
    this->D.scale((this->youngModulus * (1.0 - this->poissonRatio)) / ((1.0 + this->poissonRatio) * (1.0 - 2.0 * this->poissonRatio)));
  } else {
    throw std::runtime_error("Error in FEMdataKeeper::CreateMatrixD(float, float)");
  }
}

void FEMdataKeeper::Parse() {
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
