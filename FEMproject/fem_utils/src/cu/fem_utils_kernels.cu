#include "fem_utils_kernels.h"

#define THREADS_NUM 256

__global__
void kernelGenerateMask(float *elements, float *mask, size_t dim, size_t n) {
  size_t id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id < n) {
    size_t index = 6 * (dim - 1) * id;
    if (dim == 2) {
      mask[index + 0] = elements[3 * id + 0] * dim + 0;
      mask[index + 1] = elements[3 * id + 0] * dim + 1;
      mask[index + 2] = elements[3 * id + 1] * dim + 0;
      mask[index + 3] = elements[3 * id + 1] * dim + 1;
      mask[index + 4] = elements[3 * id + 2] * dim + 0;
      mask[index + 5] = elements[3 * id + 2] * dim + 1;
    } else if (dim == 3) {
      mask[index + 0] = elements[4 * id + 0] * dim + 0;
      mask[index + 1] = elements[4 * id + 0] * dim + 1;
      mask[index + 2] = elements[4 * id + 0] * dim + 2;
      mask[index + 3] = elements[4 * id + 1] * dim + 0;
      mask[index + 4] = elements[4 * id + 1] * dim + 1;
      mask[index + 5] = elements[4 * id + 1] * dim + 2;
      mask[index + 6] = elements[4 * id + 2] * dim + 0;
      mask[index + 7] = elements[4 * id + 2] * dim + 1;
      mask[index + 8] = elements[4 * id + 2] * dim + 2;
      mask[index + 9] = elements[4 * id + 3] * dim + 0;
      mask[index + 10] = elements[4 * id + 3] * dim + 1;
      mask[index + 11] = elements[4 * id + 3] * dim + 2;
    }
  }
}

void genMask_Ker(float *mask, float *elements, size_t dim, size_t elementsCount) {
  kernelGenerateMask<<<(elementsCount + THREADS_NUM - 1) / THREADS_NUM, THREADS_NUM>>>(elements, mask, dim, elementsCount);
}

__global__
void kernelTransformWithMask(size_t n, const float *mask, const float *src, float *dest) {
  size_t id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id < n) {
    dest[id] = src[static_cast<int>(mask[id])];
  }
}

void transformWithMask_Ker(size_t size, float *dest, const float *src, const float *mask) {
  kernelTransformWithMask<<<(size + THREADS_NUM - 1) / THREADS_NUM, THREADS_NUM>>>(size, mask, src, dest);
}

__global__
void kernelGenCoordinates(size_t n, float *coordinates, float *nodes, float *elements, size_t dim) {
  size_t id = blockIdx.x * blockDim.x + threadIdx.x;

  if (id < n) {
    size_t s = dim * (dim + 1);
    for (size_t i = 0 ; i < dim; ++i) {
      for (size_t j = 0; j < dim + 1; ++j) {
        coordinates[j + i * (dim + 1) + id * s] = nodes[i + static_cast<int>(elements[j + id * (dim + 1)]) * dim];
      }
    }
  }
}

void genCoordinates_Ker(size_t elementsCount, float *coordinates, float *nodes, float *elements, size_t dim) {
  kernelGenCoordinates<<<(elementsCount + THREADS_NUM - 1) / THREADS_NUM, THREADS_NUM>>>(elementsCount, coordinates, nodes, elements, dim);
}

__global__
void kernelGenFCoordinates(size_t n, float *fcoordinates, float *nodes, float *boundaryNodes, size_t dim) {
  size_t id = blockIdx.x * blockDim.x + threadIdx.x;

  if (id < n) {
    size_t s = dim * dim;
    float *fcoords = fcoordinates + id * s;
    for (size_t i = 0 ; i < dim; ++i) {
      for (size_t j = 0; j < dim; ++j) {
        fcoords[j + i * dim] = nodes[j + static_cast<int>(boundaryNodes[i + id * (dim + 1)]) * dim];
      }
    }
  }
}

void genFCoordinates_Ker(size_t bEdgesCount, float *fcoordinates, float *nodes, float *boundaryNodes, size_t dim) {
  kernelGenFCoordinates<<<(bEdgesCount + THREADS_NUM - 1) / THREADS_NUM, THREADS_NUM>>>(bEdgesCount, fcoordinates, nodes, boundaryNodes, dim);
}

__device__
float kernel_det3(float a0, float a1, float a2,
             float a3, float a4, float a5,
             float a6, float a7, float a8) {
  return a0 * a4 * a8 +
      a1 * a6 * a5 +
      a2 * a3 * a7 -
      a6 * a4 * a2 -
      a0 * a5 * a7 -
      a1 * a3 * a8;
}

__device__
float kernel_det3x3(const float *c) {
  return c[0] * c[4] * c[8] +
        c[1] * c[6] * c[5] +
        c[2] * c[3] * c[7] -
        c[6] * c[4] * c[2] -
        c[0] * c[5] * c[7] -
        c[1] * c[3] * c[8];
}

__device__
float kernel_det4x4(const float *m) {
  float v1 = kernel_det3(m[5], m[6], m[7], m[9],
      m[10], m[11], m[13], m[14], m[15]);
  float v2 = kernel_det3(m[1], m[2], m[3], m[9],
      m[10], m[11], m[13], m[14], m[15]);
  float v3 = kernel_det3(m[1], m[2], m[3], m[5],
      m[6], m[7], m[13], m[14], m[15]);
  float v4 = kernel_det3(m[1], m[2], m[3], m[5],
      m[6], m[7], m[9], m[10], m[11]);
  return v1 - v2 + v3 - v4;
}

__global__
void kernelCalculateArea(size_t n, const float *coords, float *Ccoords, float *areas, size_t dim) {
  size_t id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id < n) {
    float *C = Ccoords + id * (dim + 1) * (dim + 1);
    for (size_t i = 0; i < dim + 1; ++i) {
      for (size_t j = 0; j < dim + 1; ++j) {
        C[j + i * (dim + 1)] =
            j == 0 ? 1.f : coords[i + (j - 1) * (dim + 1) + id * dim * (dim + 1)];
      }
    }

    if (dim == 2) {
      areas[id] = fabs(kernel_det3x3(C));
    } else if (dim == 3) {
      areas[id] = fabs(kernel_det4x4(C));
    }
  }
}

void calculateArea_Ker(size_t size, const float *coordinates, float *Ccoords, float *areas, size_t dim) {
  kernelCalculateArea<<<(size + THREADS_NUM - 1) / THREADS_NUM, THREADS_NUM>>>(size, coordinates, Ccoords, areas, dim);
}

__global__
void kernelGenGradientMatrix(size_t n, float *Blocals, const float *coordinates, size_t dim) {
  size_t id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id < n) {
    const float *coords = coordinates + dim * (dim + 1) * id;
    float *B = Blocals + 3 * (dim - 1) * 6 * (dim - 1) * id;
    size_t col = 6 * (dim - 1);

    B[0 + col * 0] = coords[1 + (dim + 1) * 1] - coords[2 + (dim + 1) * 1];
    B[1 + col * 0] = 0.f;
    B[2 + col * 0] = coords[2 + (dim + 1) * 1] - coords[0 + (dim + 1) * 1];
    B[3 + col * 0] = 0.f;
    B[4 + col * 0] = coords[0 + (dim + 1) * 1] - coords[1 + (dim + 1) * 1];
    B[5 + col * 0] = 0.f;

    B[0 + col * 1] = 0.f;
    B[1 + col * 1] = coords[2 + (dim + 1) * 0] - coords[1 + (dim + 1) * 0];
    B[2 + col * 1] = 0.f;
    B[3 + col * 1] = coords[0 + (dim + 1) * 0] - coords[2 + (dim + 1) * 0];
    B[4 + col * 1] = 0.f;
    B[5 + col * 1] = coords[1 + (dim + 1) * 0] - coords[0 + (dim + 1) * 0];

    B[0 + col * 2] = coords[2 + (dim + 1) * 0] - coords[1 + (dim + 1) * 0];
    B[1 + col * 2] = coords[1 + (dim + 1) * 1] - coords[2 + (dim + 1) * 1];
    B[2 + col * 2] = coords[0 + (dim + 1) * 0] - coords[2 + (dim + 1) * 0];
    B[3 + col * 2] = coords[2 + (dim + 1) * 1] - coords[0 + (dim + 1) * 1];
    B[4 + col * 2] = coords[1 + (dim + 1) * 0] - coords[0 + (dim + 1) * 0];
    B[5 + col * 2] = coords[0 + (dim + 1) * 1] - coords[1 + (dim + 1) * 1];
  }
}

void genGradientMatrix_Ker(size_t elementsCount, float *Blocals, const float *coordinates, size_t dim) {
  kernelGenGradientMatrix<<<(elementsCount + THREADS_NUM - 1) / THREADS_NUM, THREADS_NUM>>>(elementsCount, Blocals, coordinates, dim);
}

__device__
size_t kernelGetId(size_t start, size_t id, size_t dim) {
  return (start + id) % (dim + 1);
}

__global__
void kernelGenGradientMatrix3D(size_t n, float *Blocals, const float *coordinates, size_t dim) {
  size_t id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id < n) {
    const float *coords = coordinates + dim * (dim + 1) * id;
    float *B = Blocals + 3 * (dim - 1) * 6 * (dim - 1) * id;
    size_t col = 6 * (dim - 1);


    for (size_t i = 0; i < dim + 1; ++i) {
      size_t start_id = i * dim;

      float b = pow(-1.f, i) * kernel_det3(1.f, coords[kernelGetId(i, 3, dim) + 1 * (dim + 1)], coords[kernelGetId(i, 3, dim) + 2 * (dim + 1)],
                                           1.f, coords[kernelGetId(i, 2, dim) + 1 * (dim + 1)], coords[kernelGetId(i, 2, dim) + 2 * (dim + 1)],
                                           1.f, coords[kernelGetId(i, 1, dim) + 1 * (dim + 1)], coords[kernelGetId(i, 1, dim) + 2 * (dim + 1)]);

      float c = pow(-1.f, i) * kernel_det3(coords[kernelGetId(i, 3, dim) + 0 * (dim + 1)], 1.f, coords[kernelGetId(i, 3, dim) + 2 * (dim + 1)],
                                           coords[kernelGetId(i, 2, dim) + 0 * (dim + 1)], 1.f, coords[kernelGetId(i, 2, dim) + 2 * (dim + 1)],
                                           coords[kernelGetId(i, 1, dim) + 0 * (dim + 1)], 1.f, coords[kernelGetId(i, 1, dim) + 2 * (dim + 1)]);

      float d = pow(-1.f, i) * kernel_det3(coords[kernelGetId(i, 3, dim) + 0 * (dim + 1)], coords[kernelGetId(i, 3, dim) + 1 * (dim + 1)], 1.f,
                                           coords[kernelGetId(i, 2, dim) + 0 * (dim + 1)], coords[kernelGetId(i, 2, dim) + 1 * (dim + 1)], 1.f,
                                           coords[kernelGetId(i, 1, dim) + 0 * (dim + 1)], coords[kernelGetId(i, 1, dim) + 1 * (dim + 1)], 1.f);

      // x
      B[start_id + 0 + 0 * col] = b;
      B[start_id + 1 + 3 * col] = b;
      B[start_id + 2 + 5 * col] = b;

      // y
      B[start_id + 0 + 3 * col] = c;
      B[start_id + 1 + 1 * col] = c;
      B[start_id + 2 + 4 * col] = c;

      // z
      B[start_id + 0 + 5 * col] = d;
      B[start_id + 1 + 4 * col] = d;
      B[start_id + 2 + 2 * col] = d;
    }
  }
}

void genGradientMatrix3D_Ker(size_t elementsCount, float *Blocals, const float *coordinates, size_t dim) {
  kernelGenGradientMatrix3D<<<(elementsCount + THREADS_NUM - 1) / THREADS_NUM, THREADS_NUM>>>(elementsCount, Blocals, coordinates, dim);
}

__global__
void kernelApplyConstraints(size_t n, const float *elements,
                            size_t constraintsCount, const float *constraints,
                            float *Klocals, size_t dim, size_t elemSize) {
  size_t id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id < n) {
    size_t s = 6 * (dim - 1);
    float *K = Klocals + s * s * id;
    const float *el = elements + elemSize * id;
    for (size_t c_id = 0; c_id < constraintsCount; ++c_id)
      for (size_t i = 0; i < elemSize; ++i)
        for (size_t j = 0; j < elemSize; ++j)
          for (size_t ilocal = 0; ilocal < dim; ++ilocal)
            for (size_t jlocal = 0; jlocal < dim; ++jlocal)
              if (dim * static_cast<int>(el[i]) + ilocal == static_cast<int>(constraints[c_id]) ||
                  dim * static_cast<int>(el[j]) + jlocal == static_cast<int>(constraints[c_id]))
                if (dim * static_cast<int>(el[i]) + ilocal != dim * static_cast<int>(el[j]) + jlocal)
                  K[s * (dim * i + ilocal) + dim * j + jlocal] = 0.f;
  }
}

void applyConstraints_Ker(size_t elementsCount, const float *elements, size_t constraintsCount, const float *constraints, float *Klocals, size_t dim, size_t elemSize) {
  kernelApplyConstraints<<<(elementsCount + THREADS_NUM - 1) / THREADS_NUM, THREADS_NUM>>>(elementsCount, elements, constraintsCount, constraints, Klocals, dim, elemSize);
}

__global__
void kernelGetDiagonalElements(size_t n, float *diag, const float *Klocals, size_t dim) {
  size_t id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id < n) {
    size_t s = 6 * (dim - 1);
    const float *K = Klocals + s * s * id;
    float *d = diag + s * id;
    for (size_t i = 0; i < s; ++i) {
      d[i] = K[i + s * i];
    }
  }
}

void getDiagonalElements_Ker(size_t elementsCount, float *diag, const float *Klocals, size_t dim) {
  kernelGetDiagonalElements<<<(elementsCount + THREADS_NUM - 1) / THREADS_NUM, THREADS_NUM>>>(elementsCount, diag, Klocals, dim);
}

__global__
void kernelCalculateLength(size_t n, const float *fcoordinates, float *bEdgesLengths, size_t dim) {
  size_t id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id < n) {
    const float *fcoords = fcoordinates + id * dim * dim;
    bEdgesLengths[id] = sqrt((fcoords[2] - fcoords[0]) * (fcoords[2] - fcoords[0]) +
                             (fcoords[3] - fcoords[1]) * (fcoords[3] - fcoords[1]));
  }
}

void calculateLength_Ker(size_t bEdgeCount, const float *fcoordinates, float *bEdgesLengths, size_t dim) {
  kernelCalculateLength<<<(bEdgeCount + THREADS_NUM - 1) / THREADS_NUM, THREADS_NUM>>>(bEdgeCount, fcoordinates, bEdgesLengths, dim);
}

__global__
void kernelCalculateLength3D(size_t n, const float *fcoordinates, float *bEdgesLengths, size_t dim) {
  size_t id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id < n) {
    const float *fcoords = fcoordinates + id * dim * dim;
    float ax = fcoords[3] - fcoords[0];
    float ay = fcoords[4] - fcoords[1];
    float az = fcoords[5] - fcoords[2];

    float bx = fcoords[6] - fcoords[0];
    float by = fcoords[7] - fcoords[1];
    float bz = fcoords[8] - fcoords[2];

    float a1 = (ay * bz - az * by);
    float a2 = (az * bx - ax * bz);
    float a3 = (ax * by - ay * bx);

    bEdgesLengths[id] = 0.5f * (a1 * a1 + a2 * a2 + a3 * a3);
  }
}

void calculateLength3D_Ker(size_t bEdgeCount, const float *fcoordinates, float *bEdgesLengths, size_t dim) {
  kernelCalculateLength3D<<<(bEdgeCount + THREADS_NUM - 1) / THREADS_NUM, THREADS_NUM>>>(bEdgeCount, fcoordinates, bEdgesLengths, dim);
}

__device__
int kernelGetLocalId(const float *elements, size_t elementsCount,
                     size_t elementId, size_t nodeId, size_t elemSize) {
  for (size_t i = 0; i < elemSize; ++i) {
    if (elements[i + elementId * elemSize] == nodeId)
      return i;
  }
  printf("Error wrong index!\n");
  return -1;
}

__global__
void kernelCalcPressures(size_t n, float *Flocals,
                    const float *boundaryAdjElems,
                    const float *boundaryNodes,
                    const float *boundaryPressures,
                    const float *bEdgesLengths,
                    const float *boundaryNormals,
                    const float *elements,
                    const float *constraintsTypes,
                    size_t elemSize, size_t dim,
                    size_t nodesCount, size_t elementsCount) {
  size_t id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id < n) {
    float coeff = dim == 2 ? -0.5f : (-1.f / 12.f);
    size_t elementId = boundaryAdjElems[id];
    for (size_t i = 0; i < dim; ++i) {
      size_t nodeId = boundaryNodes[i + id * (dim + 1)];
      int localId = kernelGetLocalId(elements, elementsCount, elementId, nodeId, elemSize);
      for (size_t j = 0; j < dim; ++j) {
        Flocals[(dim * localId + j) + elementId * 6 * (dim - 1)] = static_cast<size_t>(constraintsTypes[nodeId]) & j ?
              0.f : coeff * boundaryPressures[id] * bEdgesLengths[id] * boundaryNormals[j + dim * id];
      }
    }
  }
}

void calcPressures_Ker(size_t bEdgeCount, float *Flocals,
                         const float *boundaryAdjElems,
                         const float *boundaryNodes,
                         const float *boundaryPressures,
                         const float *bEdgesLengths,
                         const float *boundaryNormals,
                         const float *elements,
                         const float *constraintsTypes,
                         size_t elemSize, size_t dim,
                         size_t nodesCount, size_t elementsCount) {
  kernelCalcPressures<<<(bEdgeCount + THREADS_NUM - 1) / THREADS_NUM, THREADS_NUM>>>(bEdgeCount, Flocals,
                                                           boundaryAdjElems,
                                                           boundaryNodes,
                                                           boundaryPressures,
                                                           bEdgesLengths,
                                                           boundaryNormals,
                                                           elements,
                                                           constraintsTypes,
                                                           elemSize, dim,
                                                           nodesCount, elementsCount);
}

__global__
void kernelCalcLoads(size_t n, float *loads, const float *loadsNodes, const float *loadsVectors, size_t dim) {
  size_t id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id < n) {
    for (size_t d = 0; d < dim; ++d) {
      loads[d + dim * static_cast<int>(loadsNodes[id])] = loadsVectors[d + id * dim];
    }
  }
}

void calcLoads_Ker(size_t loadsCount, float *loads,
                   const float *loadsNodes, const float *loadsVectors, size_t dim) {
  kernelCalcLoads<<<(loadsCount + THREADS_NUM - 1) / THREADS_NUM, THREADS_NUM>>>(loadsCount, loads, loadsNodes, loadsVectors, dim);
}

__global__
void kernelCalculateMlocals(size_t elementsCount, bool isLumped,
                            float *Mlocals, size_t dim,
                            float rho, const float *elementsAreas,
                            const float *coordinates) {
  size_t id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id < elementsCount) {
    // TODO: add implicit scheme
    const float *coords = coordinates + dim * (dim + 1) * id;
    float *M = Mlocals + id * 6 * (dim - 1);

    float area = 0.5f * abs((coords[0] - coords[2]) *
                          (coords[4] - coords[5]) -
                          (coords[1] - coords[2]) *
                          (coords[3] - coords[5]));

    float mass = rho * area;//elementsAreas[id];

    if (isLumped) {
      M[0] = mass / 3;
      M[1] = mass / 3;
      M[2] = mass / 3;
      M[3] = mass / 3;
      M[4] = mass / 3;
      M[5] = mass / 3;
    } else {
      M[0 + 0 * 6] = mass / 6;
      M[1 + 1 * 6] = mass / 6;
      M[2 + 2 * 6] = mass / 6;
      M[3 + 3 * 6] = mass / 6;
      M[4 + 4 * 6] = mass / 6;
      M[5 + 5 * 6] = mass / 6;

      M[0 + 2 * 6] = mass / 12;
      M[1 + 3 * 6] = mass / 12;
      M[0 + 4 * 6] = mass / 12;
      M[1 + 5 * 6] = mass / 12;
      M[2 + 4 * 6] = mass / 12;
      M[3 + 5 * 6] = mass / 12;
      M[2 + 0 * 6] = mass / 12;
      M[3 + 1 * 6] = mass / 12;
      M[4 + 0 * 6] = mass / 12;
      M[5 + 1 * 6] = mass / 12;
      M[4 + 2 * 6] = mass / 12;
      M[5 + 3 * 6] = mass / 12;
    }
  }
}

void calculateMlocals_Ker(size_t elementsCount, bool isLumped, float *Mlocals, size_t dim, float rho, const float *elementsAreas, const float *coordinates) {
  kernelCalculateMlocals<<<(elementsCount + THREADS_NUM - 1) / THREADS_NUM, THREADS_NUM>>>(elementsCount, isLumped, Mlocals, dim, rho, elementsAreas, coordinates);
}

__global__
void kernelCalculateDiag(size_t elementsCount, float *diag, const float *Mlocals, const float *Klocals,
                         size_t dim, float cM, float cK, float cC, float dampAlpha, float dampBeta) {
  size_t id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id < elementsCount) {
    size_t s = 6 * (dim - 1);
    for (size_t j = 0; j < s; ++j) {
      diag[j + id * s] = cK * Klocals[j + j * s + id * s * s];
      if (cM != 0.f) {
        diag[j + id * s] += cM * Mlocals[j + j * s + id * s * s] +
                            cC * (dampAlpha * Mlocals[j + j * s + id * s * s] + dampBeta * Klocals[j + j * s + id * s * s]);

      }
    }
  }
}

void calculateDiag_Ker(size_t elementsCount, float *diag, const float *Mlocals, const float *Klocals,
                   size_t dim, float cM, float cK, float cC, float dampAlpha, float dampBeta) {
  kernelCalculateDiag<<<(elementsCount + THREADS_NUM - 1) / THREADS_NUM, THREADS_NUM>>>(elementsCount, diag, Mlocals, Klocals,
                                                            dim, cM, cK, cC, dampAlpha, dampBeta);
}
