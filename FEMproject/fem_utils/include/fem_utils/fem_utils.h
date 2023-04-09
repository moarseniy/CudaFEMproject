#ifndef FEM_UTILS_H
#define FEM_UTILS_H

#include <cassert>
#include <stdexcept>

#include <iostream>
#include <algorithm>
#include <numeric>
#include <map>
#include <set>

#include <matrix_pack/matrix_pack.h>
#include <datakeeper.h>

class ElementsData {
public:
  ElementsData(size_t DIM, size_t elementsCount, size_t nodesCount, size_t boundaryEdgesCount, DEVICE_NAME device);
  virtual ~ElementsData() {};

  virtual void getDiagonalElements() = 0;
  virtual void transformWithMask(Matrix &src, Matrix &dest) = 0;
  virtual void reductionWithMask(Matrix &src, Matrix &dest) = 0;
  virtual void reductionWithMaskAndTransform(Matrix &src, Matrix &dest, size_t size) = 0;
  virtual void applyConstraints() = 0;
  virtual void calculateKlocal() = 0;
  virtual void calculateKlocals() = 0;
  virtual void calculateFlocal() = 0;
  virtual void calculateFlocals() = 0;
  virtual void calculateLength() = 0;
  virtual void calculateArea() = 0;
  virtual void genMask() = 0;
  virtual void genAdjElements() = 0;
  virtual void genCoordinates() = 0;
  virtual void genFCoordinates() = 0;
  virtual void genGradientMatrix() = 0;

  DEVICE_NAME get_device() const;
  size_t get_dim() const;
  size_t get_elementsCount() const;

  Matrix* get_Klocals() const;
  Matrix* get_Flocals() const;
  Matrix* get_Blocals() const;
  Matrix* get_coordinates() const;
  Matrix* get_mask() const;
  Matrix* get_adjElements() const;

  Matrix* get_diag() const;
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
  DEVICE_NAME _device;

  Matrix *elements;
  Matrix *nodes;
  Matrix *constraintsIds;
  Matrix *constraintsTypes;
  Matrix *boundaryAdjElems;
  Matrix *boundaryNodes;
  Matrix *boundaryNormals;
  Matrix *boundaryPressures;
  Matrix *D;

  Matrix *Flocals; // Right part of SLAE
  Matrix *Klocals; // Stiffness Matrix
  Matrix *Blocals; // Gradient Matrix
  Matrix *tBlocals; // Transposed Gradient Matrix
  Matrix *Clocals;
  Matrix *mask;
  Matrix *mask_sorted;
  Matrix *adjElements;
  Matrix *elementsAreas;
  Matrix *bEdgesLengths;

  Matrix *coordinates;
  Matrix *fcoordinates;

  //TODO: check how to rewrite
  Matrix *diag, *r, *m, *z, *s, *p, *u, *x;
};

class CPU_ElementsData : public ElementsData {
public:
  CPU_ElementsData();
  CPU_ElementsData(const dataKeeper &dk);
  ~CPU_ElementsData();

  void getDiagonalElements() override;
  void transformWithMask(Matrix &src, Matrix &dest) override;
  void reductionWithMask(Matrix &src, Matrix &dest) override;
  void reductionWithMaskAndTransform(Matrix &src, Matrix &dest, size_t size) override;
  void applyConstraints() override;
  void calculateKlocal() override;
  void calculateKlocals() override;
  int getLocalId(size_t elementId, size_t nodeId);
  void calculateFlocal() override;
  void calculateFlocals() override;
  void calculateArea() override;
  void calculateLength() override;
  void genMask() override;
  void genAdjElements() override;
  void genCoordinates() override;
  void genFCoordinates() override;
  void genGradientMatrix() override;
};

#endif // FEM_UTILS_H
