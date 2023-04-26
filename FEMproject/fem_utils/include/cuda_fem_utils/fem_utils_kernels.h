#ifndef FEM_UTILS_KERNELS_H_
#define FEM_UTILS_KERNELS_H_

#include <cstdio>
#include <cuda_runtime.h>
#include <cuda.h>

void genMask_Ker(float *mask, float *elements, size_t dim, size_t elementsCount);
void transformWithMask_Ker(size_t size, float *dest, const float *src, const float *mask);

void genCoordinates_Ker(size_t elementsCount, float *coordinates, float *nodes, float *elements, size_t dim);
void genFCoordinates_Ker(size_t bEdgesCount, float *fcoordinates, float *nodes, float *boundaryNodes, size_t dim);

void applyConstraints_Ker(size_t elementsCount, const float *elements, size_t constraintsCount, const float *constraints, float *Klocals, size_t dim, size_t elemSize);

void calculateLength_Ker(size_t bEdgeCount, const float *fcoordinates, float *bEdgesLengths, size_t dim);
void calculateLength3D_Ker(size_t bEdgeCount, const float *fcoordinates, float *bEdgesLengths, size_t dim);
void calculateArea_Ker(size_t size, const float *coordinates, float *Ccoords, float *areas, size_t dim);

void genGradientMatrix_Ker(size_t elementsCount, float *Blocals, const float *coordinates, size_t dim);
void genGradientMatrix3D_Ker(size_t elementsCount, float *Blocals, const float *coordinates, size_t dim);
void getDiagonalElements_Ker(size_t elementsCount, float *diag, const float *Klocals, size_t dim);

void calculateFlocal_Ker(size_t bEdgeCount, float *Flocals,
                         const float *boundaryAdjElems,
                         const float *boundaryNodes,
                         const float *boundaryPressures,
                         const float *bEdgesLengths,
                         const float *boundaryNormals,
                         const float *elements,
                         const float *constraintsTypes,
                         size_t elemSize, size_t dim,
                         size_t nodesCount, size_t elementsCount);

void calculateMlocals_Ker(size_t elementsCount, bool isLumped, float *Mlocals, size_t dim, float rho, const float *elementsAreas);
void calculateDiag_Ker(size_t elementsCount, float *diag, const float *Mlocals, const float *Klocals,
                       size_t dim, float cM, float cK, float cC, float dampAlpha, float dampBeta);
#endif /* FEM_UTILS_KERNELS_H_ */
