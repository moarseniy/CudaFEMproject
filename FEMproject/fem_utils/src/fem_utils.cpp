
#include <fem_utils/fem_utils.h>
#ifdef WITH_CUDA
#include <cuda_fem_utils/cuda_fem_utils.h>
#endif

ElementsData::ElementsData(size_t DIM, size_t elementsCount, size_t boundaryEdgesCount, DEVICE_NAME device) :
  _DIM(DIM),
  _elementsCount(elementsCount),
  _boundaryEdgesCount(boundaryEdgesCount),
  _device(device) {}

ElementsData* ElementsData::setElementsData(DEVICE_NAME device, const dataKeeper &dk) {
  if (device == CPU) {
    return new CPU_ElementsData(dk);
#ifdef WITH_CUDA
  } else if (device == CUDA) {
    return new CUDA_ElementsData(dk);
#endif
  } else {
    throw std::runtime_error("ElementsData::setElementsData -> No type for device");
  }
}

DEVICE_NAME ElementsData::get_device() const {
  return _device;
}

size_t ElementsData::get_dim() const {
  return _DIM;
}

size_t ElementsData::get_elementsCount() const {
  return _elementsCount;
}

Matrix* ElementsData::get_Klocals() const {
  return Klocals;
}

Matrix* ElementsData::get_Flocals() const {
  return Flocals;
}

Matrix* ElementsData::get_Blocals() const {
  return Blocals;
}

Matrix* ElementsData::get_coordinates() const {
  return coordinates;
}

Matrix* ElementsData::get_mask() const {
  return mask;
}

Matrix* ElementsData::get_adjElements() const {
  return adjElements;
}

Matrix* ElementsData::get_diag() const {
  return diag;
}

Matrix* ElementsData::get_r() const {
  return r;
}

Matrix* ElementsData::get_m() const {
  return m;
}

Matrix* ElementsData::get_z() const {
  return z;
}

Matrix* ElementsData::get_s() const {
  return s;
}

Matrix* ElementsData::get_p() const {
  return p;
}

Matrix* ElementsData::get_u() const {
  return u;
}

Matrix* ElementsData::get_x() const {
  return x;
}
