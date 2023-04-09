#include <cuda_fem_utils/cuda_fem_utils.h>


CUDA_ElementsData::CUDA_ElementsData() :
  ElementsData(0, 0, 0, 0, CUDA) {}

CUDA_ElementsData::CUDA_ElementsData(const dataKeeper &dk) :
  ElementsData(dk.get_dim(),
               dk.get_elementsCount(),
               dk.get_nodesCount(),
               dk.get_boundaryEdgesCount(),
               CUDA) {
  nodes = Matrix::setMatrix(_device);
  dk.get_nodes()->copy(*nodes);

  elements = Matrix::setMatrix(_device);
  dk.get_elementsIds()->copy(*elements);

  constraintsIds = Matrix::setMatrix(_device);
  dk.get_constraintsIds()->copy(*constraintsIds);
  constraintsTypes = Matrix::setMatrix(_device);
  dk.get_constraintsTypes()->copy(*constraintsTypes);

  boundaryNodes = Matrix::setMatrix(_device);
  dk.get_boundaryNodes()->copy(*boundaryNodes);
  boundaryAdjElems = Matrix::setMatrix(_device);
  dk.get_boundaryAdjElems()->copy(*boundaryAdjElems);
  boundaryNormals = Matrix::setMatrix(_device);
  dk.get_boundaryNormals()->copy(*boundaryNormals);
  boundaryPressures = Matrix::setMatrix(_device);
  dk.get_boundaryPressures()->copy(*boundaryPressures);

  D = Matrix::setMatrix(_device);
  dk.get_D()->copy(*D);

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
                                  6 * (_DIM - 1));
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

CUDA_ElementsData::~CUDA_ElementsData() {
  // TODO: ADD If and check all deletess
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

  delete elementsAreas;
  delete bEdgesLengths;

  delete diag;
  delete r;
  delete m;
  delete z;
  delete s;
  delete p;
  delete u;
  delete x;
}

void CUDA_ElementsData::genMask() {
  genMask_Ker(mask->get_data(), elements->get_data(), _DIM, _elementsCount);
  mask->copy(*mask_sorted);
  mask_sorted->sort();
  genAdjElements();
}

void CUDA_ElementsData::genAdjElements() {
  assert(mask_sorted->get_numElements() > 0);
  CUDA_Matrix ones(mask_sorted->get_numRows(), mask_sorted->get_numCols());
  ones.setTo(1.f);
  ones.reduce_by_key(*mask_sorted, *adjElements);
}

void CUDA_ElementsData::getDiagonalElements() {
  getDiagonalElements_Ker(_elementsCount, diag->get_data(), Klocals->get_data(), _DIM);
}

void CUDA_ElementsData::transformWithMask(Matrix &src, Matrix &dest) {
  assert(src.get_numElements() > 0);
  transformWithMask_Ker(mask->get_numElements(), dest.get_data(), src.get_data(), mask->get_data());
}

void CUDA_ElementsData::reductionWithMask(Matrix &src, Matrix &dest) {
  Matrix *temp_mask = Matrix::setMatrix(*mask);
  mask->copy(*temp_mask);

  src.sort_by_key(*temp_mask);
  src.reduce_by_key(*temp_mask, dest);

  delete temp_mask;
}

void CUDA_ElementsData::reductionWithMaskAndTransform(Matrix &src, Matrix &dest, size_t size) {
  Matrix *temp_res = Matrix::setVector(CUDA, size);
  reductionWithMask(src, *temp_res);
  transformWithMask(*temp_res, dest);
  delete temp_res;
}

void CUDA_ElementsData::applyConstraints() {
  applyConstraints_Ker(_elementsCount, elements->get_data(), constraintsIds->get_numElements(), constraintsIds->get_data(), Klocals->get_data(), _DIM, _DIM + 1);
}

void CUDA_ElementsData::calculateArea() {
  calculateArea_Ker(_elementsCount, coordinates->get_data(), Clocals->get_data(), elementsAreas->get_data(), _DIM);
}

void CUDA_ElementsData::calculateLength() {
  calculateLength_Ker(_boundaryEdgesCount, fcoordinates->get_data(), bEdgesLengths->get_data(), _DIM);
}

void CUDA_ElementsData::genGradientMatrix() {
  genGradientMatrix_Ker(_elementsCount, Blocals->get_data(), coordinates->get_data(), _DIM);
  Blocals->divideElementwise(*elementsAreas, X);
}

void CUDA_ElementsData::genCoordinates() {
  genCoordinates_Ker(_elementsCount, coordinates->get_data(), nodes->get_data(), elements->get_data(), _DIM);
}

void CUDA_ElementsData::genFCoordinates() {
  genFCoordinates_Ker(_boundaryEdgesCount, fcoordinates->get_data(), nodes->get_data(), boundaryNodes->get_data(), _DIM);
}

void CUDA_ElementsData::calculateKlocal() {

  // K = B^T * D * B

  tBlocals->bmm(*D, 3 * (_DIM - 1), 3 * (_DIM - 1), false,
               *Blocals, 6 * (_DIM - 1), true, _elementsCount);

  Klocals->bmm(*Blocals, 6 * (_DIM - 1), 3 * (_DIM - 1), false,
               *tBlocals, 6 * (_DIM - 1), false, _elementsCount, 0.5f);

  Klocals->scale(*elementsAreas, X);
}

void CUDA_ElementsData::calculateKlocals() {
  genCoordinates();
  calculateArea();
  genGradientMatrix();

  calculateKlocal();

  applyConstraints();
  getDiagonalElements();
}

// TODO: MAKE THIS LIGHT!!!
void CUDA_ElementsData::calculateFlocal() {
  calculateFlocal_Ker(_boundaryEdgesCount, Flocals->get_data(),
                      boundaryAdjElems->get_data(),
                      boundaryNodes->get_data(),
                      boundaryPressures->get_data(),
                      bEdgesLengths->get_data(),
                      boundaryNormals->get_data(),
                      elements->get_data(),
                      constraintsTypes->get_data(),
                      _DIM + 1, _DIM,
                      _nodesCount, _elementsCount);
}

void CUDA_ElementsData::calculateFlocals() {
  genFCoordinates();
  calculateLength();
  calculateFlocal();
  CPU_Matrix tmp;
  Flocals->copy(tmp);
  tmp.Show();
  exit(-1);
}

