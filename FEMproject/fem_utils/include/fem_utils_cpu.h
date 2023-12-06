#ifndef FEM_UTILS_CPU_H
#define FEM_UTILS_CPU_H

#include <cassert>
#include <stdexcept>

#include <iostream>
#include <algorithm>
#include <numeric>
#include <map>
#include <set>

#include <matrix_pack.h>
#include <datakeeper.h>

struct CGMData {
  CGMData(DEVICE_NAME devType, size_t n_elems, size_t DIM);
  ~CGMData();

  void zeroData();
  Matrix *r, *m, *z, *s, *p, *u, *x;
};

class ElementsData {
public:
  ElementsData(size_t DIM, size_t elementsCount, size_t nodesCount, size_t boundaryEdgesCount, size_t loadsCount, DEVICE_NAME device);
  virtual ~ElementsData() {};

  float updateWavelet(float t, const WaveletParams &waveParams);
  float Ricker(float t, float ampl, float freq);
  float Berlage(float t, float ampl, float freq);

  virtual void getDiagonalElements(Matrix &Locals, Matrix &tgt) = 0;
  virtual void transformWithMask(Matrix &src, Matrix &dest) = 0;
  virtual void reductionWithMask(Matrix &src, Matrix &dest) = 0;
  virtual void reductionWithMaskAndTransform(Matrix &src, Matrix &dest, size_t size) = 0;
  virtual void applyConstraints() = 0;
  virtual void calculateKlocal() = 0;
  virtual void calculateKlocals() = 0;

  virtual void initLoads(float t, const WaveletParams &waveParams) = 0;
  virtual void calcLoads(float t, const WaveletParams &waveParams) = 0;

  virtual void initPressures(float t, const WaveletParams &waveParams) = 0;
  virtual void calcPressures(float t, const WaveletParams &waveParams) = 0;

  virtual void initFlocals(float t, const WaveletParams &waveParams) = 0;
  virtual void calcFlocals(float t, const WaveletParams &waveParams) = 0;


  virtual void calculateLength() = 0;
  virtual void calculateLength3D() = 0;
  virtual void calculateArea() = 0;
  virtual void genMask() = 0;
  virtual void genAdjElements() = 0;
  virtual void genCoordinates() = 0;
  virtual void genFCoordinates() = 0;
  virtual void genGradientMatrix() = 0;
  virtual void genGradientMatrix3D() = 0;

  virtual void calculateDiag(Matrix &diag, float cM = 0.f, float cK = 0.f, float cC = 0.f, float dampAlpha = 0.f, float dampBeta = 0.f) = 0;
  // TODO: take it to another new class
  virtual void solveDiagSystem(Matrix &diagonal, Matrix &v, Matrix &tgt, bool transformRes) = 0;
  virtual void calculateMlocals(bool isLumped, const MechanicalParams &mechParams) = 0;

  DEVICE_NAME get_device() const;
  size_t get_dim() const;
  size_t get_elementsCount() const;
  size_t get_nodesCount() const;

  Matrix* get_Klocals() const;
  Matrix* get_Flocals() const;
  Matrix* get_Blocals() const;
  Matrix* get_tBlocals() const;
  Matrix* get_D() const;
  Matrix* get_coordinates() const;
  Matrix* get_fcoordinates() const;
  Matrix* get_mask() const;
  Matrix* get_adjElements() const;
  Matrix* get_bEdgesLengths() const;
  Matrix* get_Ccoords() const;

  Matrix* get_diagK() const;

  Matrix* get_diagM() const;
  Matrix* get_Mlocals() const;
  Matrix* get_Clocals() const;

  Matrix* get_elementsAreas() const;

  void zeroCGMData();

  Matrix* get_r() const;
  Matrix* get_m() const;
  Matrix* get_z() const;
  Matrix* get_s() const;
  Matrix* get_p() const;
  Matrix* get_u() const;
  Matrix* get_x() const;

  static ElementsData* setElementsData(DEVICE_NAME device, const dataKeeper &dk);

protected:
  size_t _DIM;
  size_t _elementsCount;
  size_t _nodesCount;
  size_t _boundaryEdgesCount;
  size_t _loadsCount;
  DEVICE_NAME _device;

  Matrix *elements;
  Matrix *nodes;

  Matrix *constraintsIds;
  Matrix *constraintsTypes;

  Matrix *boundaryAdjElems;
  Matrix *boundaryNodes;
  Matrix *boundaryNormals;
  Matrix *boundaryPressures;

  Matrix *loads;
  Matrix *loadsNodes;
  Matrix *loadsVectors;

  Matrix *D;

  Matrix *Flocals; // Local Right parts of SLAE
  Matrix *Klocals; // Local Stiffness Matrices
  Matrix *diagK; // main diagonals of Stiffness Matrices
  Matrix *Blocals; // Gradient Matrices
  Matrix *tBlocals; // Temporary Gradient Matrices
  Matrix *Ccoords;
  Matrix *mask;
  Matrix *mask_sorted;
  Matrix *adjElements;
  Matrix *elementsAreas;
  Matrix *bEdgesLengths;

  Matrix *coordinates;
  Matrix *fcoordinates;

  CGMData _cgmData;

  Matrix *Mlocals;
  Matrix *diagM;
  Matrix *Clocals;
};

class CPU_ElementsData : public ElementsData {
public:
  CPU_ElementsData();
  CPU_ElementsData(const dataKeeper &dk);
  ~CPU_ElementsData();

  void getDiagonalElements(Matrix &Locals, Matrix &tgt) override;
  void transformWithMask(Matrix &src, Matrix &dest) override;
  void reductionWithMask(Matrix &src, Matrix &dest) override;
  void reductionWithMaskAndTransform(Matrix &src, Matrix &dest, size_t size) override;
  void applyConstraints() override;
  void calculateKlocal() override;
  void calculateKlocals() override;
  int getLocalId(size_t elementId, size_t nodeId);

  void initLoads(float t, const WaveletParams &waveParams) override;
  void calcLoads(float t, const WaveletParams &waveParams) override;

  void initPressures(float t, const WaveletParams &waveParams) override;
  void calcPressures(float t, const WaveletParams &waveParams) override;

  void initFlocals(float t, const WaveletParams &waveParams) override;
  void calcFlocals(float t, const WaveletParams &waveParams) override;

  void calculateArea() override;
  void calculateLength() override;
  void calculateLength3D() override;
  void genMask() override;
  void genAdjElements() override;
  void genCoordinates() override;
  void genFCoordinates() override;
  void genGradientMatrix() override;
  void genGradientMatrix2();
  void genGradientMatrix2D();
  void genGradientMatrix3D() override;

  void calculateDiag(Matrix &diag, float cM = 0.f, float cK = 0.f, float cC = 0.f, float dampAlpha = 0.f, float dampBeta = 0.f) override;

  void solveDiagSystem(Matrix &diagonal, Matrix &v, Matrix &tgt, bool transformRes) override;
  void calculateMlocals(bool isLumped, const MechanicalParams &mechParams) override;
};

#endif // FEM_UTILS_CPU_H
