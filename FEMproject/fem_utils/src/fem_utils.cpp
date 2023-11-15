
#include <fem_utils/fem_utils.h>
#ifdef WITH_CUDA
#include <cuda_fem_utils/cuda_fem_utils.h>
#endif

CGMData::CGMData(DEVICE_NAME devType, size_t n_elems, size_t DIM) {
  r = Matrix::setMatrix(devType, n_elems, 6 * (DIM - 1));
  m = Matrix::setMatrix(devType, n_elems, 6 * (DIM - 1));
  z = Matrix::setMatrix(devType, n_elems, 6 * (DIM - 1));
  s = Matrix::setMatrix(devType, n_elems, 6 * (DIM - 1));
  p = Matrix::setMatrix(devType, n_elems, 6 * (DIM - 1));
  u = Matrix::setMatrix(devType, n_elems, 6 * (DIM - 1));
  x = Matrix::setMatrix(devType, n_elems, 6 * (DIM - 1));

  zeroData();
}

CGMData::~CGMData() {
  delete r;
  delete m;
  delete z;
  delete s;
  delete p;
  delete u;
  delete x;
}

void CGMData::zeroData() {
  r->setTo(0.f); m->setTo(0.f); z->setTo(0.f);
  s->setTo(0.f); p->setTo(0.f); u->setTo(0.f);
  x->setTo(0.f);
}

ElementsData::ElementsData(size_t DIM, size_t elementsCount, size_t nodesCount,
                           size_t boundaryEdgesCount, size_t loadsCount, DEVICE_NAME device) :
  _DIM(DIM),
  _elementsCount(elementsCount),
  _nodesCount(nodesCount),
  _boundaryEdgesCount(boundaryEdgesCount),
  _loadsCount(loadsCount),
  _device(device),
  _cgmData(device, elementsCount, DIM) {}

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

float ElementsData::Ricker(float t, float ampl, float freq) {
  float tmp = M_PI * freq * t - M_PI;
  float value = ampl * (1.f - 2.f * tmp * tmp) *
                std::exp(-1.f * tmp * tmp);
  return value;
}

// Aldridge, D. F. (1990). The Berlage wavelet. GEOPHYSICS, 55(11), 1508–1511. doi:10.1190/1.1442799
// CAE-Fidesys-4.0/preprocessor/bin/help/finite_element_model/non_exodus/time_formulas.htm
float ElementsData::Berlage(float t, float ampl, float freq) {
  float w0 = 2 * M_PI * freq;
  float w1 = w0 / sqrtf(3.f);
  float value = ampl * w1 * w1 / 4.f * std::exp(-1.f * w1 * t) * (
                std::sin(w0 * t) * (1.f / (w1 * w1 * w1) + t / (w1 * w1) - t * t / w1) -
                std::cos(w0 * t) * sqrtf(3.f) * (t * t / w1 + t / (w1 * w1)) );
  return value;
}

float ElementsData::updateWavelet(float t, const WaveletParams &waveParams) {
  if (t < waveParams.timeshift) {
    return 0.f;
  } else {
    if (waveParams.waveletType == "ricker") {
      return Ricker(t - waveParams.timeshift, waveParams.ampl, waveParams.freq);
    } else if (waveParams.waveletType == "berlage") {
      return Berlage(t - waveParams.timeshift, waveParams.ampl, waveParams.freq);
    }
  }
  return 0.f; // TODO: think what to return
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

size_t ElementsData::get_nodesCount() const {
  return _nodesCount;
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

Matrix* ElementsData::get_fcoordinates() const {
  return fcoordinates;
}

Matrix* ElementsData::get_mask() const {
  return mask;
}

Matrix* ElementsData::get_adjElements() const {
  return adjElements;
}

Matrix* ElementsData::get_bEdgesLengths() const {
  return bEdgesLengths;
}

Matrix* ElementsData::get_Ccoords() const {
  return Ccoords;
}

Matrix* ElementsData::get_diagK() const {
  return diagK;
}

Matrix* ElementsData::get_diagM() const {
  return diagM;
}

Matrix* ElementsData::get_Mlocals() const {
  return Mlocals;
}

Matrix* ElementsData::get_Clocals() const {
  return Clocals;
}

Matrix* ElementsData::get_elementsAreas() const {
  return elementsAreas;
}

void ElementsData::zeroCGMData() {
  _cgmData.zeroData();
}

Matrix* ElementsData::get_r() const {
  return _cgmData.r;
}

Matrix* ElementsData::get_m() const {
  return _cgmData.m;
}

Matrix* ElementsData::get_z() const {
  return _cgmData.z;
}

Matrix* ElementsData::get_s() const {
  return _cgmData.s;
}

Matrix* ElementsData::get_p() const {
  return _cgmData.p;
}

Matrix* ElementsData::get_u() const {
  return _cgmData.u;
}

Matrix* ElementsData::get_x() const {
  return _cgmData.x;
}
