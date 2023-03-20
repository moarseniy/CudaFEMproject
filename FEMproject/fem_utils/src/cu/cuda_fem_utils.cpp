#include <cuda_fem_utils/cuda_fem_utils.h>

CUDA_ElementsData::CUDA_ElementsData() :
  ElementsData(0, 0, 0, CUDA) {}

CUDA_ElementsData::CUDA_ElementsData(const dataKeeper &dk) :
  ElementsData(dk.get_dim(),
               dk.get_elementsCount(),
               dk.get_boundaryEdgesCount(),
               CUDA) {

  nodes->setDevice(CUDA);
//  CUDA_Matrix::copy(dk.get_nodes(), *nodes);
//  nodes = Matrix::setMatrix(*dk.get_nodes());
//  dk.get_nodes().copy(*nodes);
//  elements = Matrix::setMatrix(*dk.get_elementsIds());
//  dk.get_elementsIds().copy(*elements);
//  constraints = Matrix::setMatrix(*dk.get_constraintsIds());
//  dk.get_constraintsIds().copy(*constraints);
//  D = Matrix::setMatrix(*dk.get_D());
//  dk.get_D().copy(*D);

  Blocals = Matrix::setMatrix(CPU,
                        _elementsCount,
                        3 * (_DIM - 1) * 6 * (_DIM - 1));

  Klocals = Matrix::setMatrix(CPU,
                              _elementsCount,
                              6 * (_DIM - 1) * 6 * (_DIM - 1));

  Clocals = Matrix::setMatrix(CPU,
                              _elementsCount,
                              (_DIM + 1) * (_DIM + 1));

  Flocals = Matrix::setMatrix(CPU,
                              _elementsCount,
                              6 * (_DIM - 1));
  mask = Matrix::setMatrix(CPU,
                           _elementsCount,
                           6 * (_DIM - 1));

  mask_sorted = Matrix::setMatrix(CPU,
                           _elementsCount,
                           6 * (_DIM - 1));

  adjElements = Matrix::setMatrix(CPU,
                                  _elementsCount,
                                  6 * (_DIM - 1));
  coordinates = Matrix::setMatrix(CPU,
                                  _elementsCount,
                                  6 * (_DIM - 1));
}

CUDA_ElementsData::~CUDA_ElementsData() {

}

void CUDA_ElementsData::genMask() {
  genMask_ker(mask->get_data(), elements->get_data(), _DIM, _elementsCount);
}

void CUDA_ElementsData::genAdjElements() {
  assert(false);
}

void CUDA_ElementsData::getDiagonalElements(Matrix &K, size_t el) {
  assert(false);
}

void CUDA_ElementsData::transformWithMask(Matrix &src, Matrix &dest) {
  assert(false);
}

void CUDA_ElementsData::reductionWithMask(Matrix &src, Matrix &dest) {
  assert(false);
}

void CUDA_ElementsData::reductionWithMaskAndTransform(Matrix &src, Matrix &dest, size_t size) {
  assert(false);
}

void CUDA_ElementsData::applyConstraints(size_t el) {
  assert(false);
}

float CUDA_ElementsData::calculateArea(size_t el) {
  assert(false);
  return 0.f;
}

float CUDA_ElementsData::calculateLength(size_t elementId) {
  assert(false);
  return 0.f;
}

float CUDA_ElementsData::genGradientMatrix(size_t el) {
  assert(false);
  return 0.f;
}

void CUDA_ElementsData::genCoordinates(size_t el) {
  assert(false);
}

void CUDA_ElementsData::genFCoordinates(size_t el) {
  assert(false);
}

void CUDA_ElementsData::calculateKlocal(size_t el) {
  assert(false);
}

void CUDA_ElementsData::calculateKlocals() {
  assert(false);
}

size_t CUDA_ElementsData::getLocalId(size_t elementId, size_t nodeId) {
  assert(false);
  return 0;
}

void CUDA_ElementsData::calculateFlocal(size_t el) {
  assert(false);
}

void CUDA_ElementsData::calculateFlocals() {
  assert(false);
}

