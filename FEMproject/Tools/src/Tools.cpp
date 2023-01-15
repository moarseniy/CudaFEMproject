
#include <Tools.h>

namespace fs {

void normalizeString(std::string &s) {
  if (s.length() <= 1) {
    return;
  }

  for (size_t i = s.length() - 1; i >= 1; --i) {
    if (s[i] == '/' && s[i - 1]  == '/') {
      s.erase(i, 1);
    }
  }
}

std::string joinString(const std::string &src1,
                const std::string &src2) {
  std::string str_path = src1 + "/" + src2;
  normalizeString(str_path);
  return str_path;
}

} // namespace fs
