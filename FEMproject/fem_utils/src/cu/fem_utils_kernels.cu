#include "fem_utils_kernels.h"

// nxm * kxp
__device__
void mult(float *a, float *b, float *tgt,
          size_t n, size_t m,
          size_t k, size_t p,
          bool a_tr, size_t thread_id) {
  float a1;
  if (a_tr) {
    size_t temp = n;
    n = m;
    m = temp;
  }
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < p; j++) {
      tgt[j + i * p + n * p * thread_id] = 0.f;
      for (size_t t = 0; t < n; t++) {
        a1 = a_tr ? a[i + t * m + n * m * thread_id] :
                    a[t + i * m + n * m * thread_id];
        tgt[j + i * p + n * p * thread_id] += a1 * b[j + t * p];
      }
    }
  }
}

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
  kernelGenerateMask<<<(elementsCount + 255) / 256, 256>>>(elements, mask, dim, elementsCount);
}

__global__
void kernelTransformWithMask(size_t n, const float *mask, const float *src, float *dest) {
  size_t id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id < n) {
    dest[id] = src[int(mask[id])];
  }
}

void transformWithMask_Ker(size_t size, float *dest, const float *src, const float *mask) {
  kernelTransformWithMask<<<(size + 255) / 256, 256>>>(size, mask, src, dest);
}

__global__
void kernelGenCoordinates(size_t n, float *coordinates, float *nodes, float *elements, size_t dim) {
  size_t id = blockIdx.x * blockDim.x + threadIdx.x;
  size_t s = dim * (dim + 1);
  if (id < n) {
    for (size_t i = 0 ; i < dim; ++i) {
      for (size_t j = 0; j < dim + 1; ++j) {
        coordinates[j + i * (dim + 1) + id * s] = nodes[i + static_cast<int>(elements[j + id * (dim + 1)]) * dim];
      }
    }
  }
}

void genCoordinates_Ker(size_t elementsCount, float *coordinates, float *nodes, float *elements, size_t dim) {
  kernelGenCoordinates<<<(elementsCount + 255) / 256, 256>>>(elementsCount, coordinates, nodes, elements, dim);
}

__global__
void kernelGenFCoordinates(size_t n, float *fcoordinates, float *nodes, float *boundaryNodes, size_t dim) {
  size_t id = blockIdx.x * blockDim.x + threadIdx.x;
  size_t s = dim * dim;
  if (id < n) {
    for (size_t i = 0 ; i < dim; ++i) {
      for (size_t j = 0; j < dim; ++j) {
        fcoordinates[j + i * dim + id * s] = nodes[j + static_cast<int>(boundaryNodes[i + id * dim]) * dim];
      }
    }
  }
}

void genFCoordinates_Ker(size_t bEdgesCount, float *fcoordinates, float *nodes, float *boundaryNodes, size_t dim) {
  kernelGenFCoordinates<<<(bEdgesCount + 255) / 256, 256>>>(bEdgesCount, fcoordinates, nodes, boundaryNodes, dim);
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
void kernelCalculateArea(size_t n, const float *coords, float *Clocals, float *areas, size_t dim) {
  size_t id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id < n) {
    float *C = Clocals + id * (dim + 1) * (dim + 1);
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

void calculateArea_Ker(size_t size, const float *coordinates, float *Clocals, float *areas, size_t dim) {
  kernelCalculateArea<<<(size + 255) / 256, 256>>>(size, coordinates, Clocals, areas, dim);
}

__global__
void kernelGenGradientMatrix(size_t n, float *Blocals, const float *coordinates, size_t dim) {
  size_t id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id < n) {
    const float *coords = coordinates + dim * (dim + 1) * id;
    float *B = Blocals + 3 * (dim - 1) * 6 * (dim - 1) * id;
    size_t col = 6 * (dim - 1);

    //TODO: WORKS FOR 2D ONLY!!!

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
  kernelGenGradientMatrix<<<(elementsCount + 255) / 256, 256>>>(elementsCount, Blocals, coordinates, dim);
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
  kernelApplyConstraints<<<(elementsCount + 255) / 256, 256>>>(elementsCount, elements, constraintsCount, constraints, Klocals, dim, elemSize);
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
  kernelGetDiagonalElements<<<(elementsCount + 255) / 256, 256>>>(elementsCount, diag, Klocals, dim);
}

__global__
void kernelCalculateLength(size_t n, const float *fcoordinates, float *bEdgesLengths, size_t dim) {
  size_t id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id < n) {
    const float *fcoords = fcoordinates + id * dim * dim;
    bEdgesLengths[id] = std::sqrt((fcoords[2] - fcoords[0]) * (fcoords[2] - fcoords[0]) +
                                  (fcoords[3] - fcoords[1]) * (fcoords[3] - fcoords[1]));
  }
}

void calculateLength_Ker(size_t bEdgeCount, const float *fcoordinates, float *bEdgesLengths, size_t dim) {
  kernelCalculateLength<<<(bEdgeCount + 255) / 256, 256>>>(bEdgeCount, fcoordinates, bEdgesLengths, dim);
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
void kernelCalculateFlocal(size_t n, float *Flocals,
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
    size_t elementId = boundaryAdjElems[id];
    for (size_t i = 0; i < dim; ++i) {
      size_t nodeId = boundaryNodes[i + id * (dim + 1)];
      int localId = kernelGetLocalId(elements, elementsCount, elementId, nodeId, elemSize);
      for (size_t j = 0; j < dim; ++j) {
        Flocals[(dim * localId + j) + elementId * 6 * (dim - 1)] = static_cast<size_t>(constraintsTypes[nodeId]) & j ?
              0.f : -0.5f * boundaryPressures[id] * bEdgesLengths[id] * boundaryNormals[j + dim * id];
      }
    }
  }
}

void calculateFlocal_Ker(size_t bEdgeCount, float *Flocals,
                         const float *boundaryAdjElems,
                         const float *boundaryNodes,
                         const float *boundaryPressures,
                         const float *bEdgesLengths,
                         const float *boundaryNormals,
                         const float *elements,
                         const float *constraintsTypes,
                         size_t elemSize, size_t dim,
                         size_t nodesCount, size_t elementsCount) {
  kernelCalculateFlocal<<<(bEdgeCount + 255) / 256, 256>>>(bEdgeCount, Flocals,
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
