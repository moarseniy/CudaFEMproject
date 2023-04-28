#include <cuda_fem_utils/cuda_fem_utils.h>


CUDA_ElementsData::CUDA_ElementsData() :
  ElementsData(0, 0, 0, 0, 0, CUDA) {}

CUDA_ElementsData::CUDA_ElementsData(const dataKeeper &dk) :
  ElementsData(dk.get_dim(),
               dk.get_elementsCount(),
               dk.get_nodesCount(),
               dk.get_boundaryEdgesCount(),
               dk.get_loadsCount(),
               CUDA) {
  nodes = Matrix::setMatrix(_device);
  dk.get_nodes()->copy(*nodes);

  elements = Matrix::setMatrix(_device);
  dk.get_elementsIds()->copy(*elements);

  if (dk.get_constraintsIds()->get_numElements() > 0) {
    constraintsIds = Matrix::setMatrix(_device);
    dk.get_constraintsIds()->copy(*constraintsIds);
    constraintsTypes = Matrix::setMatrix(_device);
    dk.get_constraintsTypes()->copy(*constraintsTypes);
  }

  if (_boundaryEdgesCount > 0) {
    boundaryNodes = Matrix::setMatrix(_device);
    dk.get_boundaryNodes()->copy(*boundaryNodes);
    boundaryAdjElems = Matrix::setMatrix(_device);
    dk.get_boundaryAdjElems()->copy(*boundaryAdjElems);
    boundaryNormals = Matrix::setMatrix(_device);
    dk.get_boundaryNormals()->copy(*boundaryNormals);
    boundaryPressures = Matrix::setMatrix(_device);
    dk.get_boundaryPressures()->copy(*boundaryPressures);
  }

  if (_loadsCount > 0) {
    loads = Matrix::setMatrix(_device, _nodesCount, _DIM);
    loads->setTo(0.f);
    loadsNodes = Matrix::setVector(_device);
    dk.get_loadsNodes()->copy(*loadsNodes);
    loadsVectors = Matrix::setMatrix(_device);
    dk.get_loadsVectors()->copy(*loadsVectors);
  }

  D = Matrix::setMatrix(_device);
  dk.get_D()->copy(*D);

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

CUDA_ElementsData::~CUDA_ElementsData() {
  if (nodes)
    delete nodes;
  if (elements)
    delete elements;

  if (constraintsIds)
    delete constraintsIds;
  if (constraintsTypes)
    delete constraintsTypes;

  if (_boundaryEdgesCount > 0) {
    delete boundaryAdjElems;
    delete boundaryNodes;
    delete boundaryNormals;
    delete boundaryPressures;
  }

  if (_loadsCount > 0) {
    delete loads;
    delete loadsNodes;
    delete loadsVectors;
  }

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

  if (Mlocals)
    delete Mlocals;
  if (diagM)
    delete diagM;
  if (Clocals)
    delete Clocals;
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

void CUDA_ElementsData::getDiagonalElements(Matrix &Locals, Matrix &tgt) {
  getDiagonalElements_Ker(_elementsCount, tgt.get_data(), Locals.get_data(), _DIM);
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
  calculateArea_Ker(_elementsCount, coordinates->get_data(), Ccoords->get_data(), elementsAreas->get_data(), _DIM);
}

void CUDA_ElementsData::calculateLength3D() {
  calculateLength3D_Ker(_boundaryEdgesCount, fcoordinates->get_data(), bEdgesLengths->get_data(), _DIM);
}

void CUDA_ElementsData::calculateLength() {
  calculateLength_Ker(_boundaryEdgesCount, fcoordinates->get_data(), bEdgesLengths->get_data(), _DIM);
}

void CUDA_ElementsData::genGradientMatrix3D() {
  genGradientMatrix3D_Ker(_elementsCount, Blocals->get_data(), coordinates->get_data(), _DIM);
  Blocals->divideElementwise(*elementsAreas, X);
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

  // K = B^T * D * B * Area * coeff

  float coeff = _DIM == 2 ? 0.5f : (1.f / 6.f);
//  elementsAreas->scale(coeff); // TODO: Add it to calculateArea function

  tBlocals->bmm(*D, 3 * (_DIM - 1), 3 * (_DIM - 1), false,
               *Blocals, 6 * (_DIM - 1), true, _elementsCount);

  Klocals->bmm(*Blocals, 6 * (_DIM - 1), 3 * (_DIM - 1), false,
               *tBlocals, 6 * (_DIM - 1), false, _elementsCount, coeff);

  Klocals->scale(*elementsAreas, X);
}

void CUDA_ElementsData::calculateKlocals() {
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

void CUDA_ElementsData::initFlocals(float t, const WaveletParams &waveParams) {
  if (_boundaryEdgesCount > 0) {
    initPressures(t, waveParams);
  }

//  if (_loadsCount > 0) {
//    initLoads(t, waveParams);
//  }
}

void CUDA_ElementsData::calcFlocals(float t, const WaveletParams &waveParams) {
  if (_boundaryEdgesCount > 0) {
    calcPressures(t, waveParams);
  }

//  if (_loadsCount > 0) {
//    calcLoads(t, waveParams);
//  }
}

void CUDA_ElementsData::initPressures(float t, const WaveletParams &waveParams) {
  genFCoordinates();

  if (_DIM == 3) {
    calculateLength3D();
  } else if (_DIM == 2){
    calculateLength();
  } else {
    throw std::runtime_error("wrong dimension");
  }

  calcPressures(t, waveParams);
}

// TODO: MAKE THIS LIGHT!!!
void CUDA_ElementsData::calcPressures(float t, const WaveletParams &waveParams) {
  float v = updateWavelet(t, waveParams);
  if (v != 0.f) {
    //TODO: add for each bEdge
    boundaryPressures->setTo(v);
  }

  calcPressures_Ker(_boundaryEdgesCount, Flocals->get_data(),
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

void CUDA_ElementsData::initLoads(float t, const WaveletParams &waveParams) {
  calcLoads(t, waveParams);
}

void CUDA_ElementsData::calcLoads(float t, const WaveletParams &waveParams) {
  float v = updateWavelet(t, waveParams);
  calcLoads_Ker(_loadsCount, loads->get_data(), loadsNodes->get_data(), loadsVectors->get_data(), _DIM);
  transformWithMask(*loads, *Flocals);
}

void CUDA_ElementsData::calculateMlocals(bool isLumped, const MechanicalParams &mechParams) {
  calculateMlocals_Ker(_elementsCount, isLumped, Mlocals->get_data(), _DIM, mechParams.rho, elementsAreas->get_data(), coordinates->get_data());
}

void CUDA_ElementsData::solveDiagSystem(Matrix &diagonal, Matrix &v, Matrix &tgt, bool transformRes) {
  CUDA_Matrix diagonal_assemblied(dynamic_cast<CUDA_Matrix&>(*adjElements));
  diagonal_assemblied.setTo(0.f);
  CUDA_Matrix v_assemblied(dynamic_cast<CUDA_Matrix&>(*adjElements));
  v_assemblied.setTo(0.f);

  reductionWithMask(diagonal, diagonal_assemblied);
  reductionWithMask(v, v_assemblied);

  if (transformRes) {
    v_assemblied.divideElementwise(diagonal_assemblied, tgt);
  } else {
    CUDA_Matrix temp(dynamic_cast<CUDA_Matrix&>(*adjElements));
    v_assemblied.divideElementwise(diagonal_assemblied, temp);
    transformWithMask(temp, tgt);
  }
}

void CUDA_ElementsData::calculateDiag(Matrix &diag, float cM, float cK, float cC, float dampAlpha, float dampBeta) {
  calculateDiag_Ker(_elementsCount, diag.get_data(), Mlocals->get_data(), Klocals->get_data(), _DIM, cM, cK, cC, dampAlpha, dampBeta);
}
