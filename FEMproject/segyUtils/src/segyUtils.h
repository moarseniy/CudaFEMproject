
#ifndef SEGY_UTILS_H_
#define SEGY_UTILS_H_

#include <string>
#include <vector>
#include <set>


#include <segy.h>
#include <Tools.h>

#include <matrix_pack/matrix_pack.h>
#include <cuda_matrix_pack/cuda_matrix.h>

class SegY {
public:
  SegY(int v);
  ~SegY();

  enum Type {
    Convert,
    Add
  };

  void setParameter(const std::string& filename,
    const ReceiversType::eType type, const size_t dof, Matrix &nodes,
    const size_t elementsCount, const size_t nodesCount, const size_t nt, bool isPlaneTask,
    int mpirank = 0);
  void addFilenameTXT(const std::string& txt_filename,
    const std::string& txt_filename_copy);
  void setOutputType(Type type);

  void convert() const;

private:
  std::string getLineTextHeader(const std::vector<float>& times) const;
  void fill_sgy_with_traces(segy_file* fp) const;
  void read_sgy_with_traces(segy_file* fp) const;
  int formatsize(SEGY_FORMAT format) const;

  int _mpirank;
  std::string _filename;
  ReceiversType::eType _type;
  size_t _dof;
  size_t _elementsCount;
  size_t _nodesCount;
  size_t _nt;
  bool _isPlaneTask;
  Type _output_type;

  Matrix *_nodes;

  std::set<std::string> _txt_filenames;

  static bool is_warning_shown;
};

#endif  // SEGY_UTILS_H_
