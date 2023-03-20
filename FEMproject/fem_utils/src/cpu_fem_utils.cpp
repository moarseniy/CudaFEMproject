#include <fem_utils/fem_utils.h>

CPU_ElementsData::CPU_ElementsData() :
  ElementsData(0, 0, 0, CPU) {}

CPU_ElementsData::CPU_ElementsData(const dataKeeper &dk) :
  ElementsData(dk.get_dim(),
               dk.get_elementsCount(),
               dk.get_boundaryEdgesCount(),
               CPU) {

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

//  elements = Matrix::setMatrix(*dk.get_elementsIds());
//  dk.get_elementsIds().copy(*elements);
//  constraints = Matrix::setMatrix(*dk.get_constraintsIds());
//  dk.get_constraintsIds().copy(*constraints);
//  D = Matrix::setMatrix(*dk.get_D());
//  dk.get_D().copy(*D);

  Blocals = Matrix::setMatrix(CPU,
                        _elementsCount,
                        3 * (_DIM - 1) * 6 * (_DIM - 1));

  tBlocals = Matrix::setMatrix(CPU,
                        _elementsCount,
                        6 * (_DIM - 1) * 3 * (_DIM - 1));

  Klocals = Matrix::setMatrix(CPU,
                              _elementsCount,
                              6 * (_DIM - 1) * 6 * (_DIM - 1));

  Clocals = Matrix::setMatrix(CPU,
                              _elementsCount,
                              (_DIM + 1) * (_DIM + 1));

  Flocals = Matrix::setMatrix(CPU,
                              _elementsCount,
                              6 * (_DIM - 1));
  Flocals->setTo(0.f);

  mask = Matrix::setMatrix(CPU,
                           _elementsCount,
                           6 * (_DIM - 1));

  mask_sorted = Matrix::setMatrix(CPU,
                           _elementsCount,
                           6 * (_DIM - 1));

  adjElements = Matrix::setMatrix(CPU,
                                  _DIM,
                                  dk.get_nodesCount());
  coordinates = Matrix::setMatrix(CPU,
                                  _elementsCount,
                                  6 * (_DIM - 1));
  fcoordinates = Matrix::setMatrix(CPU,
                                  dk.get_boundaryEdgesCount(),
                                  _DIM * _DIM);
  diag = Matrix::setMatrix(CPU, _elementsCount, 6 * (_DIM - 1));
  r = Matrix::setMatrix(CPU, _elementsCount, 6 * (_DIM - 1));
  m = Matrix::setMatrix(CPU, _elementsCount, 6 * (_DIM - 1));
  z = Matrix::setMatrix(CPU, _elementsCount, 6 * (_DIM - 1));
  s = Matrix::setMatrix(CPU, _elementsCount, 6 * (_DIM - 1));
  p = Matrix::setMatrix(CPU, _elementsCount, 6 * (_DIM - 1));
  u = Matrix::setMatrix(CPU, _elementsCount, 6 * (_DIM - 1));
  x = Matrix::setMatrix(CPU, _elementsCount, 6 * (_DIM - 1));
}

CPU_ElementsData::~CPU_ElementsData() {
  delete nodes;
  delete elements;
  delete constraintsIds;
  delete constraintsTypes;
  delete boundaryAdjElems;
  delete boundaryNodes;
  delete boundaryNormals;
  delete boundaryPressures;
  delete D;

  delete Blocals;
  delete tBlocals;
  delete Klocals;
  delete Clocals;
  delete Flocals;
  delete mask;
  delete mask_sorted;
  delete adjElements;

  delete coordinates;
  delete fcoordinates;

  delete diag;
  delete r;
  delete m;
  delete z;
  delete s;
  delete p;
  delete u;
  delete x;
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

void CPU_ElementsData::genCoordinates(size_t el) {
  std::unique_ptr<Matrix> coords = coordinates->getRow(el);
  coords->resize(_DIM, _DIM + 1);
  for (size_t i = 0 ; i < _DIM; ++i) {
    for (size_t j = 0; j < _DIM + 1; ++j) {
      (*coords)(i, j) = (*nodes)(int((*elements)(el, j)), i);
    }
  }
//  for (size_t i = 0 ; i < _DIM + 1; ++i) {
//    for (size_t j = 0; j < _DIM; ++j) {
//      (*coords)[k] = (*nodes)(int((*elements)(el, i)), j);
//      k++;
//    }
//  }
}

// should be depricated
void CPU_ElementsData::genFCoordinates(size_t el) {
  std::unique_ptr<Matrix> coords = fcoordinates->getRow(el);
  coords->resize(_DIM, _DIM);
  for (size_t i = 0 ; i < _DIM; ++i) {
    for (size_t j = 0; j < _DIM; ++j) {
      (*coords)(i, j) = (*nodes)(int((*boundaryNodes)(el, i)), j);
    }
  }
}

float CPU_ElementsData::calculateLength(size_t el) {
  std::unique_ptr<Matrix> coords = fcoordinates->getRow(el);
  return std::sqrt(((*coords)[2] - (*coords)[0]) * ((*coords)[2] - (*coords)[0]) + ((*coords)[3] - (*coords)[1]) * ((*coords)[3] - (*coords)[1]));
}

float CPU_ElementsData::calculateArea(size_t el) {
  std::unique_ptr<Matrix> C = Clocals->getRow(el);
  std::unique_ptr<Matrix> coords = coordinates->getRow(el);
//  std::cout << coords->get_numElements() << "\n";
  coords->resize(_DIM, _DIM + 1);
  C->resize(_DIM + 1, _DIM + 1);
//  coords->resize(_DIM + 1, _DIM + 1);
  for (size_t i = 0; i < _DIM + 1; ++i) {
    for (size_t j = 0; j < _DIM + 1; ++j) {
      C->get_data()[j + i * (_DIM + 1)] =
          j == 0 ? 1.f : (*coords)(j - 1, i);//[(j - 1) + i * (_DIM)];
    }
  }
//  C->Show();
  return C->det();
}

float CPU_ElementsData::genGradientMatrix(size_t el) {
  float area = calculateArea(el);
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

  B->divideElementwise(std::abs(area));
  return area;
}

void CPU_ElementsData::reductionWithMask(Matrix &src, Matrix &dest) {
  Matrix *temp_mask = Matrix::setMatrix(*mask);
  mask->copy(*temp_mask);

  src.sort_by_key(*temp_mask);
  src.reduce_by_key(*temp_mask, dest);

  delete temp_mask;
}

void CPU_ElementsData::transformWithMask(Matrix &src, Matrix &dest) {
  for (size_t i = 0; i < dest.get_numElements(); ++i) {
    dest[i] = src[int((*mask)[i])];
  }
}

void CPU_ElementsData::reductionWithMaskAndTransform(Matrix &src, Matrix &dest, size_t size) {
  Matrix *temp_res = Matrix::setVector(CPU, size);
  reductionWithMask(src, *temp_res);
  transformWithMask(*temp_res, dest);
  delete temp_res;
}

void CPU_ElementsData::applyConstraints(size_t el) {
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

void CPU_ElementsData::getDiagonalElements(Matrix &K, size_t el) {
  std::unique_ptr<Matrix> diagonal = diag->getRow(el);
  for (size_t i = 0; i < K.get_numRows(); ++i) {
    (*diagonal)[i] = K(i, i);
  }
}

void CPU_ElementsData::calculateKlocal(size_t el) {
  genCoordinates(el);
  float area = genGradientMatrix(el);

  std::unique_ptr<Matrix> B = Blocals->getRow(el);
  std::unique_ptr<Matrix> temp = tBlocals->getRow(el);
  std::unique_ptr<Matrix> K = Klocals->getRow(el);
  B->resize(3 * (_DIM - 1), 6 * (_DIM - 1));
  temp->resize(6 * (_DIM - 1), 3 * (_DIM - 1));
  K->resize(6 * (_DIM - 1), 6 * (_DIM - 1));

  B->product(*D, *temp, true);
  temp->product(*B, *K);
  K->scale(std::abs(area) * 0.5f);

  getDiagonalElements(*K, el);
}

void CPU_ElementsData::calculateKlocals() {
  for (size_t el = 0; el < _elementsCount; ++el) {
    calculateKlocal(el);
    applyConstraints(el);
  }
//  Klocals->Show();
}

size_t CPU_ElementsData::getLocalId(size_t elementId, size_t nodeId) {
  for (size_t i = 0; i < _DIM + 1; ++i) {
    if ((*elements)(elementId, i) == nodeId)
      return i;
  }
  return -1;
}

void CPU_ElementsData::calculateFlocal(size_t bEdge) {
  size_t elementId = (*boundaryAdjElems)[bEdge];
  genFCoordinates(bEdge);
  float edge_length = calculateLength(bEdge);
//  std::cout << "len = " << edge_length << "\n";
  for (size_t i = 0; i < _DIM; ++i) {
    size_t nodeId = (*boundaryNodes)(bEdge, i);
    size_t localId = getLocalId(elementId, nodeId);
    for (size_t j = 0; j < _DIM; ++j) {
      (*Flocals)(elementId, _DIM * localId + j) = int((*constraintsTypes)(nodeId, i)) & static_cast<Constraint::Type>(j) ?
            0.f : -0.5f * (*boundaryPressures)[bEdge] * edge_length * (*boundaryNormals)(bEdge, j);
    }
  }
}

void CPU_ElementsData::calculateFlocals() {
  for (size_t bEdge = 0; bEdge < _boundaryEdgesCount; ++bEdge) {
    calculateFlocal(bEdge);
  }
//  Flocals->Show();
}
