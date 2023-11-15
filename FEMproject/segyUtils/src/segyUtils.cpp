#include "segyUtils.h"

#include <locale>
#include <codecvt>

#include <fstream>
#include <sstream>
#include <cstdio>
#include <chrono>
#include <ctime>

std::string getLineTextHeader(const std::vector<float>& times, const ReceiversType::eType type, const size_t dof, size_t elementsCount, size_t nodesCount) {
  std::string line_textheader =
    "C 1 PROGRAM: CudaFEMproject                                                       "
    "C 2 LICENSE:                                                                    "
    "C 3 CLIENT:                                                                     "
    "C 4 RECORDING PARAMETER:                                                        "
    "C 5 FILE FORMAT: SEG-Y Rev.1                                                    "
    "C 6 DATE:                                                                       "
    "C 7 NUMBER OF ELEMENTS IN CALCULATED MODEL:                                     "
    "C 8 NUMBER OF NODES IN CALCULATED MODEL:                                        "
    "C 9 NUMBER OF SAMPLES IN TRACE:                                                 "
    "C10                                                                             "
    "C11 FIRST SAMPLE:                                                               "
    "C12 LAST SAMPLE:                                                                "
    "C13 SAMPLE INTERVAL:                                                            "
    "C14                                                                             "
    "C15 TRACE HEADER POSITION:                                                      "
    "C16 OFFSET BYTES 037-041                                                        "
    "C17 INLINE BYTES 189-193                                                        "
    "C18 CROSSLINE BYTES 193-197                                                     "
    "C19                                                                             "
    "C20 END EBCDIC HEADER                                                           "
    "C21                                                                             "
    "C22                                                                             "
    "C23                                                                             "
    "C24                                                                             "
    "C25                                                                             "
    "C26                                                                             "
    "C27                                                                             "
    "C28                                                                             "
    "C29                                                                             "
    "C30                                                                             "
    "C31                                                                             "
    "C32                                                                             "
    "C33                                                                             "
    "C34                                                                             "
    "C35                                                                             "
    "C36                                                                             "
    "C37                                                                             "
    "C38                                                                             "
    "C39                                                                             "
    "C40                                                                            \x80";
  const std::string version_str = "1.0";
  strcpy(&line_textheader[0 * 80 + 25], version_str.c_str());

  const std::string license_str = "Unprotected";
  const std::string client_str = "CudaFEMproject";

  strcpy(&line_textheader[1 * 80 + 13], license_str.c_str());
  strcpy(&line_textheader[2 * 80 + 12], client_str.c_str());
  const std::string type_str = ReceiversTypeToInfoName(type, dof);
  strcpy(&line_textheader[3 * 80 + 25], type_str.c_str());
  const auto current = std::chrono::system_clock::now();
  const auto time = std::chrono::system_clock::to_time_t(current);
  std::string time_str = std::ctime(&time);
  time_str.resize(time_str.size() - 1);
  strcpy(&line_textheader[5 * 80 + 10], time_str.c_str());
  const std::string elements_str = std::to_string(elementsCount);
  strcpy(&line_textheader[6 * 80 + 44], elements_str.c_str());
  const std::string nodes_str = std::to_string(nodesCount);
  strcpy(&line_textheader[7 * 80 + 41], nodes_str.c_str());
  const std::string timesteps_str = std::to_string(times.size());
  strcpy(&line_textheader[8 * 80 + 33], timesteps_str.c_str());
  if (times.size() > 0) {
    int first_sample = static_cast<int>(times[0] * 1e3);
    first_sample += (times[0] * 1e3 - first_sample) > 0 ? 1 : 0;
    const std::string first_sample_str = std::to_string(first_sample) + " MS";
    strcpy(&line_textheader[10 * 80 + 18], first_sample_str.c_str());
    int last_sample = static_cast<int>(times[times.size() - 1] * 1e3);
    last_sample += (times[times.size() - 1] * 1e3 - last_sample) > 0 ? 1 : 0;
    const std::string last_sample_str = std::to_string(last_sample) + " MS";
    strcpy(&line_textheader[11 * 80 + 17], last_sample_str.c_str());
    strcpy(&line_textheader[12 * 80 + 21], first_sample_str.c_str());
  }

  return line_textheader;
}

void fill_sgy_with_traces(segy_file* fp, std::string txt_filename, const ReceiversType::eType type, Matrix &nodes, size_t dof, size_t elementsCount, size_t nodesCount) {
  size_t total_trace_num(0);
  size_t line_counter(1);
  bool is_warning_shown = false;
//  for (const auto& txt_filename : _txt_filenames) {
    std::ifstream in(txt_filename, std::ifstream::in);
    if (in.good()) {
      std::string str;
      std::getline(in, str);
      std::istringstream idss(str);
      size_t trace_num;
      idss >> trace_num;

      std::vector<size_t> trace_ids(trace_num);
      size_t id_counter(0);
      size_t num;
      while (idss >> num) {
        trace_ids[id_counter] = num;
        id_counter++;
      }

      std::vector<std::vector<float>> traces(trace_num);
      std::vector<float> times;
      size_t sample_counter(0);
      while (std::getline(in, str)) {
        std::istringstream ss(str);
        float val;
        ss >> val;
        times.push_back(val);
        size_t trace_counter(0);
        while (ss >> val) {
          traces[trace_counter].push_back(static_cast<float>(val));
          trace_counter++;
        }
        sample_counter++;
      }

      const int format = SEGY_IBM_FLOAT_4_BYTE;

      int err(0);

      const std::string line_textheader = getLineTextHeader(times, type, dof, elementsCount, nodesCount);
      segy_write_textheader(fp, 0, line_textheader.c_str());

      std::string line_binheader;
      line_binheader.resize(SEGY_BINARY_HEADER_SIZE);
      //char line_binheader[SEGY_BINARY_HEADER_SIZE];
      err = segy_set_bfield(&line_binheader[0], SEGY_BINFIELD::SEGY_BIN_TRACES, static_cast<int>(total_trace_num + trace_num));
      size_t interval(0);
      if (times.size() > 0) {
        interval = static_cast<int>(times[0] * 1e6);
        interval += (times[0] * 1.0e6 - static_cast<float>(interval)) > 0 ? 1 : 0;
      }
      if (!is_warning_shown && interval > std::numeric_limits<short int>::max()) {
        std::cerr << "Timestep is too big to output in sgy-file. The max value is " + std::to_string(std::numeric_limits<short int>::max()) + " microseconds, sample time is converted to milliseconds.\n";
        interval /= 1000;
        is_warning_shown = true;
      }
      err = segy_set_bfield(&line_binheader[0], SEGY_BINFIELD::SEGY_BIN_INTERVAL, static_cast<int>(interval));
      err = segy_set_bfield(&line_binheader[0], SEGY_BINFIELD::SEGY_BIN_SAMPLES, static_cast<int>(sample_counter));
      err = segy_set_bfield(&line_binheader[0], SEGY_BINFIELD::SEGY_BIN_FORMAT, format);
      err = segy_set_bfield(&line_binheader[0], SEGY_BINFIELD::SEGY_BIN_SEGY_REVISION, 1);
      err = segy_set_bfield(&line_binheader[0], SEGY_BINFIELD::SEGY_BIN_SORTING_CODE, 4);
      err = segy_set_bfield(&line_binheader[0], SEGY_BINFIELD::SEGY_BIN_MEASUREMENT_SYSTEM, 1);
      err = segy_set_bfield(&line_binheader[0], SEGY_BINFIELD::SEGY_BIN_TRACE_FLAG, 1);
      err = segy_write_binheader(fp, line_binheader.c_str());

      for (size_t i = 0; i < static_cast<size_t>(traces.size()); ++i) {
        const int trace_bsize = static_cast<int>(traces[i].size() * sizeof(float));
        const int traceno = static_cast<int>(total_trace_num + i);
        const int trace0 = SEGY_TEXT_HEADER_SIZE + SEGY_BINARY_HEADER_SIZE;
        std::string trace_binheader;
        trace_binheader.resize(SEGY_TRACE_HEADER_SIZE);
        const size_t trace_id = trace_ids[i];
        if (trace_id < 0)
          throw std::runtime_error("ERROR: Invalid trace_id in SegY::fill_sgy_with_traces.");
        err = segy_set_field(&trace_binheader[0], SEGY_FIELD::SEGY_TR_SOURCE_X, static_cast<int>(nodes(trace_id, 0)));
        err = segy_set_field(&trace_binheader[0], SEGY_FIELD::SEGY_TR_GROUP_X, static_cast<int>(nodes(trace_id, 0)));
        err = segy_set_field(&trace_binheader[0], SEGY_FIELD::SEGY_TR_CDP_X, static_cast<int>(nodes(trace_id, 0)));
        bool isPlaneTask = true;
        if (isPlaneTask) {
          err = segy_set_field(&trace_binheader[0], SEGY_FIELD::SEGY_TR_RECV_GROUP_ELEV, static_cast<int>(nodes(trace_id, 1)));
          err = segy_set_field(&trace_binheader[0], SEGY_FIELD::SEGY_TR_SOURCE_SURF_ELEV, static_cast<int>(nodes(trace_id, 1)));
        } else {
          err = segy_set_field(&trace_binheader[0], SEGY_FIELD::SEGY_TR_SOURCE_Y, static_cast<int>(nodes(trace_id, 1)));
          err = segy_set_field(&trace_binheader[0], SEGY_FIELD::SEGY_TR_GROUP_Y, static_cast<int>(nodes(trace_id, 1)));
          err = segy_set_field(&trace_binheader[0], SEGY_FIELD::SEGY_TR_CDP_Y, static_cast<int>(nodes(trace_id, 1)));
          err = segy_set_field(&trace_binheader[0], SEGY_FIELD::SEGY_TR_RECV_GROUP_ELEV, static_cast<int>(nodes(trace_id, 2)));
          err = segy_set_field(&trace_binheader[0], SEGY_FIELD::SEGY_TR_SOURCE_SURF_ELEV, static_cast<int>(nodes(trace_id, 2)));
        }
        err = segy_set_field(&trace_binheader[0], SEGY_FIELD::SEGY_TR_COORD_UNITS, 1);
        err = segy_set_field(&trace_binheader[0], SEGY_FIELD::SEGY_TR_SAMPLE_COUNT, static_cast<int>(sample_counter));
        err = segy_set_field(&trace_binheader[0], SEGY_FIELD::SEGY_TR_SAMPLE_INTER, static_cast<int>(interval));
        err = segy_set_field(&trace_binheader[0], SEGY_FIELD::SEGY_TR_INLINE, static_cast<int>(trace_id)); // _ins->nodes.id.get(trace_id)
        err = segy_set_field(&trace_binheader[0], SEGY_FIELD::SEGY_TR_CROSSLINE, static_cast<int>(line_counter));
        err = segy_set_field(&trace_binheader[0], SEGY_FIELD::SEGY_TR_SOURCE_TYPE, static_cast<int>(dof + 1));

        err = segy_write_traceheader(fp, traceno, trace_binheader.c_str(), trace0, trace_bsize);
        err = segy_set_format(fp, format | SEGY_MSB);
        err = segy_from_native(format, traces[i].size(), traces[i].data());
        err = segy_writetrace(fp, traceno, traces[i].data(), trace0, trace_bsize);
      }

      total_trace_num += trace_num;
      line_counter++;
    }
//  }
}

void convert(std::string out_path, std::string txt_filename, const ReceiversType::eType type, Matrix &nodes, size_t dof, size_t elementsCount, size_t nodesCount) {

//#ifdef _WIN32
//    using convert_type = std::codecvt_utf8<wchar_t>;
//    std::wstring_convert<convert_type, wchar_t> converter;
//    const std::string wfilename = converter.to_bytes(_filename);
//#else
  const std::string wfilename = out_path;
//#endif  // _WIN32

  segy_file* fp = segy_open(wfilename.c_str(), "wb");
  fill_sgy_with_traces(fp, txt_filename, type, nodes, dof, elementsCount, nodesCount);
  if (fp) {
    const int err = segy_close(fp);
    if (err != 0)
      throw std::runtime_error("Can't access to sgy-file.");
  }
}

