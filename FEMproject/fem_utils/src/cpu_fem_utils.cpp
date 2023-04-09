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

//  elements = Matrix::setMatrix(*dk.get_elementsIds());
//  dk.get_elementsIds().copy(*elements);
//  constraints = Matrix::setMatrix(*dk.get_constraintsIds());
//  dk.get_constraintsIds().copy(*constraints);
//  D = Matrix::setMatrix(*dk.get_D());
//  dk.get_D().copy(*D);

  Blocals = Matrix::setMatrix(_device,
                        _elementsCount,
                        3 * (_DIM - 1) * 6 * (_DIM - 1));

  tBlocals = Matrix::setMatrix(_device,
                        _elementsCount,
                        6 * (_DIM - 1) * 3 * (_DIM - 1));

  Klocals = Matrix::setMatrix(_device,
                              _elementsCount,
                              6 * (_DIM - 1) * 6 * (_DIM - 1));

  Clocals = Matrix::setMatrix(_device,
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

  diag = Matrix::setMatrix(_device, _elementsCount, 6 * (_DIM - 1));
  r = Matrix::setMatrix(_device, _elementsCount, 6 * (_DIM - 1));
  m = Matrix::setMatrix(_device, _elementsCount, 6 * (_DIM - 1));
  z = Matrix::setMatrix(_device, _elementsCount, 6 * (_DIM - 1));
  s = Matrix::setMatrix(_device, _elementsCount, 6 * (_DIM - 1));
  p = Matrix::setMatrix(_device, _elementsCount, 6 * (_DIM - 1));
  u = Matrix::setMatrix(_device, _elementsCount, 6 * (_DIM - 1));
  x = Matrix::setMatrix(_device, _elementsCount, 6 * (_DIM - 1));
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
  delete elementsAreas;
  delete bEdgesLengths;

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
//  for (size_t i = 0 ; i < _DIM + 1; ++i) {
//    for (size_t j = 0; j < _DIM; ++j) {
//      (*coords)[k] = (*nodes)(int((*elements)(el, i)), j);
//      k++;
//    }
//  }
}

// should be depricated
void CPU_ElementsData::genFCoordinates() {
  for (size_t bEdge = 0; bEdge < _boundaryEdgesCount; ++bEdge) {
    std::unique_ptr<Matrix> coords = fcoordinates->getRow(bEdge);
    coords->resize(_DIM, _DIM);
    for (size_t i = 0 ; i < _DIM; ++i) {
      for (size_t j = 0; j < _DIM; ++j) {
        (*coords)(i, j) = (*nodes)(int((*boundaryNodes)(bEdge, i)), j);
      }
    }
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
    std::unique_ptr<Matrix> C = Clocals->getRow(el);
    std::unique_ptr<Matrix> coords = coordinates->getRow(el);
    coords->resize(_DIM, _DIM + 1);
    C->resize(_DIM + 1, _DIM + 1);
    for (size_t i = 0; i < _DIM + 1; ++i) {
      for (size_t j = 0; j < _DIM + 1; ++j) {
        C->get_data()[j + i * (_DIM + 1)] =
            j == 0 ? 1.f : (*coords)(j - 1, i);//[(j - 1) + i * (_DIM)];
      }
    }
  //  C->Show();
    (*elementsAreas)[el] = std::abs(C->det());
  }
}

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

    B->divideElementwise(std::abs((*elementsAreas)[el]));
  }
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

void CPU_ElementsData::getDiagonalElements() {
  for (size_t el = 0; el < _elementsCount; ++el) {
    std::unique_ptr<Matrix> diagonal = diag->getRow(el);
    std::unique_ptr<Matrix> K = Klocals->getRow(el);
    for (size_t i = 0; i < K->get_numRows(); ++i) {
      (*diagonal)[i] = (*K)(i, i);
    }
  }
}

void CPU_ElementsData::calculateKlocal() {
  genCoordinates();
  calculateArea();
  genGradientMatrix();

  // K = B^T * D * B

  tBlocals->bmm(*D, 3 * (_DIM - 1), 3 * (_DIM - 1), false,
               *Blocals, 6 * (_DIM - 1), true, _elementsCount);

  Klocals->bmm(*Blocals, 6 * (_DIM - 1), 3 * (_DIM - 1), false,
               *tBlocals, 6 * (_DIM - 1), false, _elementsCount, 0.5f);

  Klocals->scale(*elementsAreas, X);
  getDiagonalElements();

//  Klocals->Show();

//  std::unique_ptr<Matrix> B = Blocals->getRow(el);
//  std::unique_ptr<Matrix> temp = tBlocals->getRow(el);
//  std::unique_ptr<Matrix> K = Klocals->getRow(el);
//  B->resize(3 * (_DIM - 1), 6 * (_DIM - 1));
//  temp->resize(6 * (_DIM - 1), 3 * (_DIM - 1));
//  K->resize(6 * (_DIM - 1), 6 * (_DIM - 1));

//  B->product(*D, *temp, true);
//  temp->product(*B, *K);
//  K->scale(std::abs((*elementsAreas)[el]) * 0.5f);
}

void CPU_ElementsData::calculateKlocals() {
  calculateKlocal();
  applyConstraints();
//  Klocals->Show();
//  exit(-1);
}

int CPU_ElementsData::getLocalId(size_t elementId, size_t nodeId) {
  for (size_t i = 0; i < _DIM + 1; ++i) {
    if ((*elements)(elementId, i) == nodeId)
      return i;
  }
  return -1;
}

void CPU_ElementsData::calculateFlocal() {
  genFCoordinates();
  calculateLength();
  boundaryAdjElems->Show();
  for (size_t bEdge = 0; bEdge < _boundaryEdgesCount; ++bEdge) {
    size_t elementId = (*boundaryAdjElems)[bEdge];
//    printf("%d:%d ", bEdge, elementId);
  //  std::cout << "len = " << edge_length << "\n";
    for (size_t i = 0; i < _DIM; ++i) {
      size_t nodeId = (*boundaryNodes)(bEdge, i);
      int localId = getLocalId(elementId, nodeId);
      std::cout << localId << "\n";
      for (size_t j = 0; j < _DIM; ++j) {
        (*Flocals)(elementId, _DIM * localId + j) = int((*constraintsTypes)[nodeId]) & static_cast<Constraint::Type>(j) ?
              0.f : -0.5f * (*boundaryPressures)[bEdge] * (*bEdgesLengths)[bEdge] * (*boundaryNormals)(bEdge, j);
      }
    }
  }
}

void CPU_ElementsData::calculateFlocals() {
  calculateFlocal();
  Flocals->Show();
  exit(-1);
}
