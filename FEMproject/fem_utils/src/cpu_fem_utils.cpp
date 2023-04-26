#include <fem_utils/fem_utils.h>

CPU_ElementsData::CPU_ElementsData() :
  ElementsData(0, 0, 0, 0, CPU) {}

CPU_ElementsData::CPU_ElementsData(const dataKeeper &dk) :
  ElementsData(dk.get_dim(),
               dk.get_elementsCount(),
               dk.get_nodesCount(),
               dk.get_boundaryEdgesCount(),
               CPU) {

//  dk.get_nodes()->copy(*nodes);
  nodes = Matrix::setMatrix(*dk.get_nodes());
  CPU_Matrix::copy(*dk.get_nodes(), *nodes);

  elements = Matrix::setMatrix(*dk.get_elementsIds());
  CPU_Matrix::copy(*dk.get_elementsIds(), *elements);

  constraintsIds = Matrix::setMatrix(*dk.get_constraintsIds());
  CPU_Matrix::copy(*dk.get_constraintsIds(), *constraintsIds);
  constraintsTypes = Matrix::setMatrix(*dk.get_constraintsTypes());
  CPU_Matrix::copy(*dk.get_constraintsTypes(), *constraintsTypes);

  boundaryNodes = Matrix::setMatrix(*dk.get_boundaryNodes());
  CPU_Matrix::copy(*dk.get_boundaryNodes(), *boundaryNodes);
  boundaryAdjElems = Matrix::setMatrix(*dk.get_boundaryAdjElems());
  CPU_Matrix::copy(*dk.get_boundaryAdjElems(), *boundaryAdjElems);
  boundaryNormals = Matrix::setMatrix(*dk.get_boundaryNormals());
  CPU_Matrix::copy(*dk.get_boundaryNormals(), *boundaryNormals);
  boundaryPressures = Matrix::setMatrix(*dk.get_boundaryPressures());
  CPU_Matrix::copy(*dk.get_boundaryPressures(), *boundaryPressures);

  D = Matrix::setMatrix(*dk.get_D());
  CPU_Matrix::copy(*dk.get_D(), *D);

  Blocals = Matrix::setMatrix(_device,
                        _elementsCount,
                        3 * (_DIM - 1) * 6 * (_DIM - 1));
  Blocals->setTo(0.f);

  tBlocals = Matrix::setMatrix(_device,
                        _elementsCount,
                        6 * (_DIM - 1) * 3 * (_DIM - 1));

  Klocals = Matrix::setMatrix(_device,
                              _elementsCount,
                              6 * (_DIM - 1) * 6 * (_DIM - 1));
  Klocals->setTo(0.f);

  Ccoords = Matrix::setMatrix(_device,
                              _elementsCount,
                              (_DIM + 1) * (_DIM + 1));

  Flocals = Matrix::setMatrix(_device,
                              _elementsCount,
                              6 * (_DIM - 1));
  Flocals->setTo(0.f);

  mask = Matrix::setMatrix(_device,
                           _elementsCount,
                           6 * (_DIM - 1));

  mask_sorted = Matrix::setMatrix(_device,
                           _elementsCount,
                           6 * (_DIM - 1));

  adjElements = Matrix::setMatrix(_device,
                                  _DIM,
                                  dk.get_nodesCount());
  coordinates = Matrix::setMatrix(_device,
                                  _elementsCount,
                                  _DIM * (_DIM + 1));
  fcoordinates = Matrix::setMatrix(_device,
                                  dk.get_boundaryEdgesCount(),
                                  _DIM * _DIM);

  elementsAreas = Matrix::setVector(_device, dk.get_elementsCount());
  bEdgesLengths = Matrix::setVector(_device, dk.get_boundaryEdgesCount());

  diagK = Matrix::setMatrix(_device, _elementsCount, 6 * (_DIM - 1));

//  if (dk.getTaskType() == "dynamic") {
    Mlocals = Matrix::setMatrix(_device,
                                _elementsCount,
                                6 * (_DIM - 1));

    diagM = Matrix::setMatrix(_device, _elementsCount, 6 * (_DIM - 1));

    Clocals = Matrix::setMatrix(_device,
                                _elementsCount,
                                6 * (_DIM - 1) * 6 * (_DIM - 1));

//  }
}

CPU_ElementsData::~CPU_ElementsData() {
  if (nodes)
    delete nodes;
  if (elements)
    delete elements;
  if (constraintsIds)
    delete constraintsIds;
  if (constraintsTypes)
    delete constraintsTypes;
  if (boundaryAdjElems)
    delete boundaryAdjElems;
  if (boundaryNodes)
    delete boundaryNodes;
  if (boundaryNormals)
    delete boundaryNormals;
  if (boundaryPressures)
    delete boundaryPressures;
  if (D)
    delete D;

  if (Blocals)
    delete Blocals;
  if (tBlocals)
    delete tBlocals;
  if (Klocals)
    delete Klocals;
  if (Ccoords)
    delete Ccoords;
  if (Flocals)
    delete Flocals;
  if (mask)
    delete mask;
  if (mask_sorted)
    delete mask_sorted;
  if (adjElements)
    delete adjElements;

  if (coordinates)
    delete coordinates;
  if (fcoordinates)
    delete fcoordinates;

  if (elementsAreas)
    delete elementsAreas;
  if (bEdgesLengths)
    delete bEdgesLengths;

  if (diagK)
    delete diagK;

  if (Mlocals != nullptr)
    delete Mlocals;
  if (diagM != nullptr)
    delete diagM;
  if (Clocals != nullptr)
    delete Clocals;
}

void CPU_ElementsData::genMask() {
  float *mask_ptr = mask->get_data();
  float *elements_ptr = elements->get_data();
  for (size_t i = 0; i < _elementsCount; ++i) {
    size_t index = 6 * (_DIM - 1) * i;
    if (_DIM == 2) {
      mask_ptr[index + 0] = elements_ptr[3 * i + 0] * _DIM + 0;
      mask_ptr[index + 1] = elements_ptr[3 * i + 0] * _DIM + 1;
      mask_ptr[index + 2] = elements_ptr[3 * i + 1] * _DIM + 0;
      mask_ptr[index + 3] = elements_ptr[3 * i + 1] * _DIM + 1;
      mask_ptr[index + 4] = elements_ptr[3 * i + 2] * _DIM + 0;
      mask_ptr[index + 5] = elements_ptr[3 * i + 2] * _DIM + 1;
    } else if (_DIM == 3) {
      mask_ptr[index + 0] = elements_ptr[4 * i + 0] * _DIM + 0;
      mask_ptr[index + 1] = elements_ptr[4 * i + 0] * _DIM + 1;
      mask_ptr[index + 2] = elements_ptr[4 * i + 0] * _DIM + 2;
      mask_ptr[index + 3] = elements_ptr[4 * i + 1] * _DIM + 0;
      mask_ptr[index + 4] = elements_ptr[4 * i + 1] * _DIM + 1;
      mask_ptr[index + 5] = elements_ptr[4 * i + 1] * _DIM + 2;
      mask_ptr[index + 6] = elements_ptr[4 * i + 2] * _DIM + 0;
      mask_ptr[index + 7] = elements_ptr[4 * i + 2] * _DIM + 1;
      mask_ptr[index + 8] = elements_ptr[4 * i + 2] * _DIM + 2;
      mask_ptr[index + 9] = elements_ptr[4 * i + 3] * _DIM + 0;
      mask_ptr[index + 10] = elements_ptr[4 * i + 3] * _DIM + 1;
      mask_ptr[index + 11] = elements_ptr[4 * i + 3] * _DIM + 2;
    }
  }

  mask->sort(*mask_sorted);

  genAdjElements();
}

void CPU_ElementsData::genAdjElements() {
  CPU_Matrix ones(mask_sorted->get_numRows(), mask_sorted->get_numCols());
  ones.setTo(1.f);
  ones.reduce_by_key(*mask_sorted, *adjElements);
}

void CPU_ElementsData::genCoordinates() {
  for (size_t el = 0; el < _elementsCount; ++el) {
    std::unique_ptr<Matrix> coords = coordinates->getRow(el);
    coords->resize(_DIM, _DIM + 1);
    for (size_t i = 0 ; i < _DIM; ++i) {
      for (size_t j = 0; j < _DIM + 1; ++j) {
        (*coords)(i, j) = (*nodes)(int((*elements)(el, j)), i);
      }
    }
  }
}

// should be depricated
void CPU_ElementsData::genFCoordinates() {
  for (size_t bEdge = 0; bEdge < _boundaryEdgesCount; ++bEdge) {
    std::unique_ptr<Matrix> coords = fcoordinates->getRow(bEdge);
    coords->resize(_DIM, _DIM);
    for (size_t i = 0 ; i < _DIM; ++i) {
      for (size_t j = 0; j < _DIM; ++j) {
        (*coords)(i, j) = (*nodes)(static_cast<int>((*boundaryNodes)(bEdge, i)), j);
//        std::cout << (*nodes)(static_cast<int>((*boundaryNodes)(bEdge, i)), j) << "\n";
      }
    }
  }
}

void CPU_ElementsData::calculateLength3D() {
  for (size_t bEdge = 0; bEdge < _boundaryEdgesCount; ++bEdge) {
    std::unique_ptr<Matrix> coords = fcoordinates->getRow(bEdge);
    coords->resize(_DIM, _DIM); // x0, y0, z0, x1, y1, z1, x2, y2, z2
    float ax = (*coords)[3] - (*coords)[0];
    float ay = (*coords)[4] - (*coords)[1];
    float az = (*coords)[5] - (*coords)[2];

    float bx = (*coords)[6] - (*coords)[0];
    float by = (*coords)[7] - (*coords)[1];
    float bz = (*coords)[8] - (*coords)[2];

    float a1 = (ay * bz - az * by);
    float a2 = (az * bx - ax * bz);
    float a3 = (ax * by - ay * bx);

    (*bEdgesLengths)[bEdge] = 0.5f * (a1 * a1 + a2 * a2 + a3 * a3);
  }
}

void CPU_ElementsData::calculateLength() {
  for (size_t bEdge = 0; bEdge < _boundaryEdgesCount; ++bEdge) {
    std::unique_ptr<Matrix> coords = fcoordinates->getRow(bEdge);
    (*bEdgesLengths)[bEdge] = std::sqrt(((*coords)[2] - (*coords)[0]) * ((*coords)[2] - (*coords)[0]) + ((*coords)[3] - (*coords)[1]) * ((*coords)[3] - (*coords)[1]));
  }
}

void CPU_ElementsData::calculateArea() {
  for (size_t el = 0; el < _elementsCount; ++el) {
    std::unique_ptr<Matrix> C = Ccoords->getRow(el);
    std::unique_ptr<Matrix> coords = coordinates->getRow(el);
    coords->resize(_DIM, _DIM + 1);
    C->resize(_DIM + 1, _DIM + 1);
    for (size_t i = 0; i < _DIM + 1; ++i) {
      for (size_t j = 0; j < _DIM + 1; ++j) {
        C->get_data()[j + i * (_DIM + 1)] =
            j == 0 ? 1.f : (*coords)(j - 1, i);//[(j - 1) + i * (_DIM)];
      }
    }
    (*elementsAreas)[el] = std::abs(C->det());
  }
}

float det3_cpu(float a0, float a1, float a2,
             float a3, float a4, float a5,
             float a6, float a7, float a8) {
  return a0 * a4 * a8 +
      a1 * a6 * a5 +
      a2 * a3 * a7 -
      a6 * a4 * a2 -
      a0 * a5 * a7 -
      a1 * a3 * a8;
}

void CPU_ElementsData::genGradientMatrix2D() {
  size_t el = 0;
//  for (size_t el = 0; el < _elementsCount; ++el) {
  //  std::cout << "Element #" << el << ": area = " << area << "\n";

    std::unique_ptr<Matrix> coords = coordinates->getRow(el);
    std::unique_ptr<Matrix> B = Blocals->getRow(el);
    coords->resize(_DIM, _DIM + 1);
    B->resize(3 * (_DIM - 1), 6 * (_DIM - 1));
    std::vector<std::string> coo(6);
    coo[0] = "x1"; coo[1] = "x2"; coo[2] = "x3";
    coo[3] = "y1"; coo[4] = "y2"; coo[5] = "y3";


    for (size_t i = 0; i < _DIM + 1; ++i) {
      size_t start_id = i * _DIM;

      // x
      std::cout << "B(" << start_id << "," << 1 <<  ") = ";
      std::cout << coo[(i + 2) % (_DIM + 1) + 0 * (_DIM + 1)] << " - " << coo[(i + 1) % (_DIM + 1) + 0 * (_DIM + 1)] << "\n";
      std::cout << "B(" << start_id + 1 << "," << 2 <<  ") = ";
      std::cout << coo[(i + 2) % (_DIM + 1) + 0 * (_DIM + 1)] << " - " << coo[(i + 1) % (_DIM + 1) + 0 * (_DIM + 1)] << "\n";

      // y
      std::cout << "B(" << start_id << "," << 0 << ") = ";
      std::cout << coo[(i + 1) % (_DIM + 1) + 1 * (_DIM + 1)] << " - " << coo[(i + 2) % (_DIM + 1) + 1 * (_DIM + 1)] << "\n";
      std::cout << "B(" << start_id + 1 << "," << 2 << ") = ";
      std::cout << coo[(i + 1) % (_DIM + 1) + 1 * (_DIM + 1)] << " - " << coo[(i + 2) % (_DIM + 1) + 1 * (_DIM + 1)] << "\n";
      std::cout << "\n";
    }
//  }
}

size_t get_id(size_t start, size_t id, size_t dim) {
  return (start + id) % (dim + 1);
}

void CPU_ElementsData::genGradientMatrix3D() {
  for (size_t el = 0; el < _elementsCount; ++el) {

    std::unique_ptr<Matrix> coords = coordinates->getRow(el);
    std::unique_ptr<Matrix> B = Blocals->getRow(el);
    coords->resize(_DIM, _DIM + 1);
    B->resize(3 * (_DIM - 1), 6 * (_DIM - 1));

//    std::vector<std::string> coo(12);
//    coo[0] = "x1"; coo[1] = "x2"; coo[2] = "x3"; coo[3] = "x4";
//    coo[4] = "y1"; coo[5] = "y2"; coo[6] = "y3"; coo[7] = "y4";
//    coo[8] = "z1"; coo[9] = "z2"; coo[10] = "z3"; coo[11] = "z4";

    for (size_t i = 0; i < _DIM + 1; ++i) {
      size_t start_id = i * _DIM;

//      std::cout << (i + 1) % (_DIM + 1) << "\n";
      float b = std::pow(-1.f, i) * det3_cpu(1.f, (*coords)[get_id(i, 3, _DIM) + 1 * (_DIM + 1)], (*coords)[get_id(i, 3, _DIM) + 2 * (_DIM + 1)],
                                             1.f, (*coords)[get_id(i, 2, _DIM) + 1 * (_DIM + 1)], (*coords)[get_id(i, 2, _DIM) + 2 * (_DIM + 1)],
                                             1.f, (*coords)[get_id(i, 1, _DIM) + 1 * (_DIM + 1)], (*coords)[get_id(i, 1, _DIM) + 2 * (_DIM + 1)]);

      float c = std::pow(-1.f, i) * det3_cpu((*coords)[get_id(i, 3, _DIM) + 0 * (_DIM + 1)], 1.f, (*coords)[get_id(i, 3, _DIM) + 2 * (_DIM + 1)],
                                             (*coords)[get_id(i, 2, _DIM) + 0 * (_DIM + 1)], 1.f, (*coords)[get_id(i, 2, _DIM) + 2 * (_DIM + 1)],
                                             (*coords)[get_id(i, 1, _DIM) + 0 * (_DIM + 1)], 1.f, (*coords)[get_id(i, 1, _DIM) + 2 * (_DIM + 1)]);

      float d = std::pow(-1.f, i) * det3_cpu((*coords)[get_id(i, 3, _DIM) + 0 * (_DIM + 1)], (*coords)[get_id(i, 3, _DIM) + 1 * (_DIM + 1)], 1.f,
                                             (*coords)[get_id(i, 2, _DIM) + 0 * (_DIM + 1)], (*coords)[get_id(i, 2, _DIM) + 1 * (_DIM + 1)], 1.f,
                                             (*coords)[get_id(i, 1, _DIM) + 0 * (_DIM + 1)], (*coords)[get_id(i, 1, _DIM) + 1 * (_DIM + 1)], 1.f);

      // x
//      std::cout << "B(" << start_id + 0 << "," << 0 << ") = b_" << i << "\n";
//      std::cout << "B(" << start_id + 1 << "," << 3 << ") = b_" << i << "\n";
//      std::cout << "B(" << start_id + 2 << "," << 5 << ") = b_" << i << "\n";
      (*B)(0, start_id + 0) = b;
      (*B)(3, start_id + 1) = b;
      (*B)(5, start_id + 2) = b;

      // y
//      std::cout << "B(" << start_id + 0 << "," << 3 << ") = c_" << i << "\n";
//      std::cout << "B(" << start_id + 1 << "," << 1 << ") = c_" << i << "\n";
//      std::cout << "B(" << start_id + 2 << "," << 4 << ") = c_" << i << "\n";
      (*B)(3, start_id + 0) = c;
      (*B)(1, start_id + 1) = c;
      (*B)(4, start_id + 2) = c;

      // z
//      std::cout << "B(" << start_id + 0 << "," << 5 << ") = d_" << i << "\n";
//      std::cout << "B(" << start_id + 1 << "," << 4 << ") = d_" << i << "\n";
//      std::cout << "B(" << start_id + 2 << "," << 2 << ") = d_" << i << "\n";
      (*B)(5, start_id + 0) = d;
      (*B)(4, start_id + 1) = d;
      (*B)(2, start_id + 2) = d;

//      std::cout << "\n";
//      std::cout << d / (*elementsAreas)[el] << " ";
    }
//    std::cout << "\n";
  }
  Blocals->divideElementwise(*elementsAreas, X);
}

// TODO: Calculate it using determinants based formulas
void CPU_ElementsData::genGradientMatrix() {
  for (size_t el = 0; el < _elementsCount; ++el) {
  //  std::cout << "Element #" << el << ": area = " << area << "\n";

    std::unique_ptr<Matrix> coords = coordinates->getRow(el);
    std::unique_ptr<Matrix> B = Blocals->getRow(el);
    coords->resize(_DIM, _DIM + 1);
    B->resize(3 * (_DIM - 1), 6 * (_DIM - 1));

    //TODO: WORKS FOR 2D ONLY!!!

    (*B)(0, 0) = (*coords)(1, 1) - (*coords)(1, 2);
    (*B)(0, 1) = 0.f;
    (*B)(0, 2) = (*coords)(1, 2) - (*coords)(1, 0);
    (*B)(0, 3) = 0.f;
    (*B)(0, 4) = (*coords)(1, 0) - (*coords)(1, 1);
    (*B)(0, 5) = 0.f;

    (*B)(1, 0) = 0.f;
    (*B)(1, 1) = (*coords)(0, 2) - (*coords)(0, 1);
    (*B)(1, 2) = 0.f;
    (*B)(1, 3) = (*coords)(0, 0) - (*coords)(0, 2);
    (*B)(1, 4) = 0.f;
    (*B)(1, 5) = (*coords)(0, 1) - (*coords)(0, 0);

    (*B)(2, 0) = (*coords)(0, 2) - (*coords)(0, 1);
    (*B)(2, 1) = (*coords)(1, 1) - (*coords)(1, 2);
    (*B)(2, 2) = (*coords)(0, 0) - (*coords)(0, 2);
    (*B)(2, 3) = (*coords)(1, 2) - (*coords)(1, 0);
    (*B)(2, 4) = (*coords)(0, 1) - (*coords)(0, 0);
    (*B)(2, 5) = (*coords)(1, 0) - (*coords)(1, 1);

    //B->divideElementwise((*elementsAreas)[el]); // TODO: divide all matrix after this function
  }
  Blocals->divideElementwise(*elementsAreas, X);
}

void CPU_ElementsData::reductionWithMask(Matrix &src, Matrix &dest) {
  Matrix *temp_mask = Matrix::setMatrix(*mask);
  mask->copy(*temp_mask);

  src.sort_by_key(*temp_mask);
  src.reduce_by_key(*temp_mask, dest);

  delete temp_mask;
}

void CPU_ElementsData::transformWithMask(Matrix &src, Matrix &dest) {
  assert(src.get_numElements() > 0);
  for (size_t i = 0; i < mask->get_numElements(); ++i) {
    dest[i] = src[int((*mask)[i])];
  }
}

void CPU_ElementsData::reductionWithMaskAndTransform(Matrix &src, Matrix &dest, size_t size) {
  Matrix *temp_res = Matrix::setVector(CPU, size);
  reductionWithMask(src, *temp_res);
  transformWithMask(*temp_res, dest);
  delete temp_res;
}

void CPU_ElementsData::applyConstraints() {
  for (size_t el = 0; el < _elementsCount; ++el) {
    std::unique_ptr<Matrix> Klocal = Klocals->getRow(el);
    Klocal->resize(6 * (_DIM - 1), 6 * (_DIM - 1));
    for (size_t c_id = 0; c_id < constraintsIds->get_numElements(); ++c_id)
      for (size_t i = 0; i < _DIM + 1; ++i)
        for (size_t j = 0; j < _DIM + 1; ++j)
          for (size_t ilocal = 0; ilocal < _DIM; ++ilocal)
            for (size_t jlocal = 0; jlocal < _DIM; ++jlocal)
              if (_DIM * int((*elements)(el, i)) + ilocal == int((*constraintsIds)[c_id]) ||
                  _DIM * int((*elements)(el, j)) + jlocal == int((*constraintsIds)[c_id]))
                if (_DIM * int((*elements)(el, i)) + ilocal != _DIM * int((*elements)(el, j)) + jlocal)
                  (*Klocal)(_DIM * i + ilocal, _DIM * j + jlocal) = 0.f;

  }
}

void CPU_ElementsData::getDiagonalElements(Matrix &Locals, Matrix &tgt) {
  for (size_t el = 0; el < _elementsCount; ++el) {
    std::unique_ptr<Matrix> diagonal = tgt.getRow(el);
    std::unique_ptr<Matrix> loc = Locals.getRow(el);
    loc->resize(6 * (_DIM - 1), 6 * (_DIM - 1));
    for (size_t i = 0; i < loc->get_numRows(); ++i) {
      (*diagonal)[i] = (*loc)(i, i);
    }
  }
}

void CPU_ElementsData::calculateKlocal() {

  // K = B^T * D * B * Area * Coeff

  float coeff = _DIM == 2 ? 0.5 : (1.f / 6.f);

  tBlocals->bmm(*D, 3 * (_DIM - 1), 3 * (_DIM - 1), false,
               *Blocals, 6 * (_DIM - 1), true, _elementsCount);

  Klocals->bmm(*Blocals, 6 * (_DIM - 1), 3 * (_DIM - 1), false,
               *tBlocals, 6 * (_DIM - 1), false, _elementsCount, coeff);

  Klocals->scale(*elementsAreas, X);
}

void CPU_ElementsData::calculateKlocals() {
  genCoordinates();
  calculateArea();

  if (_DIM == 3) {
    genGradientMatrix3D();
  } else if (_DIM == 2) {
    genGradientMatrix();
  } else {
    throw std::runtime_error("wrong dimension");
  }

  calculateKlocal();

  applyConstraints();
}

int CPU_ElementsData::getLocalId(size_t elementId, size_t nodeId) {
  for (size_t i = 0; i < _DIM + 1; ++i) {
    if ((*elements)(elementId, i) == nodeId)
      return i;
  }
  return -1;
}

void CPU_ElementsData::calculateFlocal(float t, const WaveletParams &waveParams) {
  updateWavelet(t, waveParams);
  float coeff = _DIM == 2 ? -0.5f : (-1.f / 12.f);
  for (size_t bEdge = 0; bEdge < _boundaryEdgesCount; ++bEdge) {
    size_t elementId = (*boundaryAdjElems)[bEdge];
    for (size_t i = 0; i < _DIM; ++i) {
      size_t nodeId = (*boundaryNodes)(bEdge, i);
      int localId = getLocalId(elementId, nodeId);
      for (size_t j = 0; j < _DIM; ++j) {
        (*Flocals)(elementId, _DIM * localId + j) = static_cast<int>((*constraintsTypes)[nodeId]) & static_cast<Constraint::Type>(j) ?
              0.f : coeff * (*boundaryPressures)[bEdge] * (*bEdgesLengths)[bEdge] * (*boundaryNormals)(bEdge, j);
      }
    }
  }
}

void CPU_ElementsData::calculateFlocals(float t, const WaveletParams &waveParams) {
  genFCoordinates();
  if (_DIM == 3) {
    calculateLength3D();
  } else if (_DIM == 2){
    calculateLength();
  } else {
    throw std::runtime_error("wrong dimension");
  }

  calculateFlocal(t, waveParams);
}

void CPU_ElementsData::calculateMlocals(bool isLumped, const MechanicalParams &mechParams) {
  // TODO: add implicit scheme
  for (size_t el = 0; el < _elementsCount; ++el) {
    std::unique_ptr<Matrix> M = Mlocals->getRow(el);
//    M->resize(6 * (_DIM - 1), 6 * (_DIM - 1));

    std::unique_ptr<Matrix> coords = coordinates->getRow(el);
    coords->resize(_DIM, _DIM + 1);
    float area = std::abs(((*coords)[0] - (*coords)[2]) *
                          ((*coords)[4] - (*coords)[5]) -
                          ((*coords)[1] - (*coords)[2]) *
                          ((*coords)[3] - (*coords)[5]));
//    if (area != (*elementsAreas)[el])
//      std::cout << "AAAAAAAAAAA\n";
//    std::cout << area << " " << (*elementsAreas)[el] << "\n";
    float mass = mechParams.rho * (*elementsAreas)[el];

    if (isLumped) {
      (*M)[0] = mass / 3;
      (*M)[1] = mass / 3;
      (*M)[2] = mass / 3;
      (*M)[3] = mass / 3;
      (*M)[4] = mass / 3;
      (*M)[5] = mass / 3;
    } else {
      (*M)(0, 0) = mass / 6;
      (*M)(1, 1) = mass / 6;
      (*M)(2, 2) = mass / 6;
      (*M)(3, 3) = mass / 6;
      (*M)(4, 4) = mass / 6;
      (*M)(5, 5) = mass / 6;

      (*M)(2, 0) = mass / 12;
      (*M)(3, 1) = mass / 12;
      (*M)(4, 0) = mass / 12;
      (*M)(5, 1) = mass / 12;
      (*M)(4, 2) = mass / 12;
      (*M)(5, 3) = mass / 12;
      (*M)(0, 2) = mass / 12;
      (*M)(1, 3) = mass / 12;
      (*M)(0, 4) = mass / 12;
      (*M)(1, 5) = mass / 12;
      (*M)(2, 4) = mass / 12;
      (*M)(3, 5) = mass / 12;
    }
  }
}

void CPU_ElementsData::solveDiagSystem(Matrix &diagonal, Matrix &v, Matrix &tgt, bool transformRes) {
  CPU_Matrix diagonal_assemblied(dynamic_cast<CPU_Matrix&>(*adjElements));
  CPU_Matrix v_assemblied(dynamic_cast<CPU_Matrix&>(*adjElements));

  reductionWithMask(diagonal, diagonal_assemblied);
  reductionWithMask(v, v_assemblied);

  if (transformRes) {
    v_assemblied.divideElementwise(diagonal_assemblied, tgt);
  } else {
    CPU_Matrix temp(dynamic_cast<CPU_Matrix&>(*adjElements));
    v_assemblied.divideElementwise(diagonal_assemblied, temp);
    transformWithMask(temp, tgt);
  }
}

void CPU_ElementsData::calculateDiag(Matrix &diag, float cM, float cK, float cC, float dampAlpha, float dampBeta) {
  for (size_t el = 0; el < _elementsCount; ++el) {
    std::unique_ptr<Matrix> diagonal = diag.getRow(el);
    std::unique_ptr<Matrix> K = Klocals->getRow(el);
    K->resize(6 * (_DIM - 1), 6 * (_DIM - 1));

    std::unique_ptr<Matrix> M;
    if (cM != 0.f) {
      M = Mlocals->getRow(el);
//      M->resize(6 * (_DIM - 1), 6 * (_DIM - 1));
    }

    for (size_t i = 0; i < K->get_numRows(); ++i) {
      (*diagonal)[i] = cK * (*K)(i, i);
      if (cM != 0.f) {
        (*diagonal)[i] += cM * (*M)[i] + cC * (dampAlpha * (*M)[i] + dampBeta * (*K)(i, i));
      }
    }
  }
}
