
#include <VTKfile.h>

void testSEGY(std::string path1, std::string path2) {

  CSegyIn input(path1);

  if (input.open()) {
      std::cout << "SUCCESS OPEN\n";
  }
  if (input.read()) {
      std::cout << "SUCCESS READ\n";
  }

  int format = input.getSegyFormat();
  int id = input.getSegyID();
  float interval = input.getSampleInterval();
  int trNum = input.getTraceNumber();
  int sampNum = input.getSampleNumber();
  float *data = new float[trNum * sampNum];
  input.getSegyTrace(data);
//  std::cout << data[0] << " " << data[1] << " " << data[2] << "\n";

  std::cout << "FORMAT: " << format << "\nID: " << id << "\nSampleInterval: " << interval << "\nTraceNumber: " << trNum << "\nSampleNumber: " << sampNum << "\n";

  CSegyOut output(path2, input.getTraceNumber());
  if (output.open())
      std::cout << "SUCCESS OPEN\n";

  output.setSegyTrace(data, trNum * sampNum);
  output.setSampleNumber(input.getSampleNumber());
  output.setTraceNumber(input.getTraceNumber());
  output.setSegyFormat(5); // TODO: add IBM_32 format
  output.setSegyID(0);
  output.setSampleInterval(input.getSampleInterval());

  if (output.write())
      std::cout << "SUCCESS WRITE\n";
  output.close();
}

void getSEGYinfo(std::string path) {
  CSegyIn in(path);
  if (!in.open()) {
    std::cout << "Problem with open\n";
  }
  if (!in.read()) {
    std::cout << "Problem with read\n";
  }

  int format = in.getSegyFormat();
  int id = in.getSegyID();
  float interval = in.getSampleInterval();
  int trNum = in.getTraceNumber();
  int sampNum = in.getSampleNumber();

  float *data = new float[trNum * sampNum];
  in.getSegyTrace(data);

  for (size_t trace_n = 0; trace_n < trNum; ++trace_n) {
    for (size_t samp_n = 0; samp_n < sampNum; ++samp_n) {
      std::cout << data[samp_n + trace_n * sampNum] << " ";
    }
    std::cout << "\n";
  }
  delete [] data;

  std::cout << "FORMAT: " << format << "\nID: " << id << "\nSampleInterval: " << interval << "\nTraceNumber: " << trNum << "\nSampleNumber: " << sampNum << "\n";
  in.close();
}

void writeDisplForSEGY(std::string out_path, Matrix &src, std::vector<int> receivers, int axis, float timestep) {
  std::fstream output;
  output.open(out_path, std::ios::app);
  output << timestep << " ";
  for (size_t i = 0; i < receivers.size(); ++i) {
    int n = receivers[i];
    std::unique_ptr<Matrix> v = src.getRow(n);
    if (axis == 0)
      output << (*v)[0] << " "; // x
    else if (axis == 1)
      output << (*v)[1] << " "; // y
  }
  output << "\n";

  output.close();
}

void convertToSEGY(std::string src_path, std::string SEGY_path, std::vector<int> receivers, float sampleInterval, int sampleNumber) {

  int traceNumber = receivers.size();

  float *segy_data = new float[traceNumber * sampleNumber];

  std::fstream f;
  f.open(src_path, std::ios::in);
  if (f) {
    for (size_t samp_n = 0; samp_n < sampleNumber; ++samp_n) {
      for (size_t trace_n = 0; trace_n < traceNumber; ++trace_n) {
        f >> segy_data[samp_n + trace_n * sampleNumber];
//        segy_data[samp_n + trace_n * sampleNumber] = 1.f;
//        std::cout << segy_data[samp_n + trace_n * sampleNumber] << " ";
      }
//      std::cout << "\n";
    }
  } else {
    std::cout << "FILE ERROR\n";
  }

  for (size_t trace_n = 0; trace_n < traceNumber; ++trace_n) {
    for (size_t samp_n = 0; samp_n < sampleNumber; ++samp_n) {
      std::cout << segy_data[samp_n + trace_n * sampleNumber] << " ";
    }
    std::cout << "\n";
  }

  CSegyOut out(SEGY_path, traceNumber);

  if (!out.open())
    std::cout << "Problem with open\n";

  if (!out.setSegyTrace(segy_data, traceNumber * sampleNumber))
    std::cout << "Problem with setSegyTrace\n";
  out.setSegyFormat(5); // fidesys: 1
  out.setSegyID(0); // fidesys: 0
  out.setTraceNumber(traceNumber);
  out.setSampleInterval(sampleInterval);
  out.setSampleNumber(sampleNumber);

  if (!out.write())
    std::cout << "Problem with write\n";
  else
    std::cout << "SEG-Y converted!\n";

  out.close();

  f.close();
  delete [] segy_data;
}

void MakeVTKfile2D(std::string output_vtk,
                   std::vector<CPU_Matrix> &nodes,
                   std::vector<Element> &elements,
                   Matrix &displacements,
                   std::vector<CPU_Matrix> &Stress,
                   std::vector<float> &sigma_mises,
                   std::vector<CPU_Matrix> &Deformation,
                   std::vector<float> &epsilon_mises,
                   CPU_Matrix &SmoothStress) {
  CheckRunTime(__func__)
  std::fstream outvtk;
  outvtk.open(output_vtk, std::fstream::out);
  outvtk << "# vtk DataFile Version 1.0\nresults.vtk  2D Unstructured Grid of Linear Triangles\nASCII\n\nDATASET UNSTRUCTURED_GRID\nPOINTS "
           << nodes[0].get_numElements() << " double\n";
  for (size_t i = 0; i < nodes[0].get_numElements(); ++i) {
    outvtk << nodes[0][i] << " " << nodes[1][i] << " " << 0.0 << "\n";
  }

  outvtk << "CELLS " << elements.size() << " " << elements.size() * 4 << "\n";
  for (size_t i = 0; i < elements.size(); ++i) {
    outvtk << 3 << " " << elements[i].nodesIds[0] << " " << elements[i].nodesIds[1] << " " << elements[i].nodesIds[2] << "\n";
  }

  outvtk << "CELL_TYPES " << elements.size() << "\n";
  for (size_t i = 0; i < elements.size(); ++i) {
    outvtk << 5 << "\n";
  }

  outvtk << "\nPOINT_DATA " << nodes[0].get_numElements() << "\n";
  outvtk << "VECTORS displacements double\n";
  for (size_t i = 0; i < displacements.get_numElements() - 1; i += 2) {
    outvtk << displacements[i] << " " << displacements[i + 1] << " 0.0\n";
  }

  outvtk << "\nSCALARS summary double\nLOOKUP_TABLE default\n";
  for (size_t i = 0; i < displacements.get_numElements() - 1; i += 2) {
    outvtk << std::sqrt(displacements[i] * displacements[i] + displacements[i + 1] * displacements[i + 1]) << "\n";
  }

//  outvtk << "\nSCALARS SmoothStress double\nLOOKUP_TABLE default\n";
//  for (size_t i = 0; i < SmoothStress.get_numElements(); ++i) {
//    outvtk << SmoothStress[i] << "\n";
//  }

//  outvtk << "\nCELL_DATA " << elements.size() << "\n";
//  outvtk << "TENSORS stress double\n";
//  for (size_t i = 0; i < Stress.size(); i++) {
//    outvtk << Stress[i][0] << " " << Stress[i][2] << " " << 0.0 << "\n";
//    outvtk << Stress[i][2] << " " << Stress[i][1] << " " << 0.0 << "\n";
//    outvtk << 0.0 << " " << 0.0 << " " << 0.0 << "\n";
//  }

//  outvtk << "\nSCALARS XY_STRESS double\nLOOKUP_TABLE default\n";
//  for (size_t i = 0; i < elements.size(); i++) {
//    outvtk <<Stress[i][2] << "\n";
//  }

//  outvtk << "\nSCALARS XY_DEFORM double\nLOOKUP_TABLE default\n";
//  for (size_t i = 0; i < elements.size(); i++) {
//    outvtk <<Deformation[i][2] << "\n";
//  }

//  outvtk << "\nSCALARS mises_stress double\nLOOKUP_TABLE default\n";
//  for (size_t i = 0; i < sigma_mises.size(); i++) {
//    outvtk << sigma_mises[i] << "\n";
//  }

//  outvtk << "\nVECTORS deformation double\n";
//  for (size_t i = 0; i < Deformation.size(); i++) {
//    outvtk << Deformation[i][0] << " " << Deformation[i][1] << " " << Deformation[i][2] << "\n";
//  }

//  outvtk << "\nSCALARS mises_deformation double\nLOOKUP_TABLE default\n";
//  for (size_t i = 0; i < epsilon_mises.size(); i++) {
//    outvtk << epsilon_mises[i] << "\n";
//  }
}

void MakeVTKfile2D(std::string output_vtk,
                   Matrix &nodes,
                   Matrix &elements,
                   Matrix &displacements) {
  CheckRunTime(__func__)
  std::fstream outvtk;
  outvtk.open(output_vtk, std::fstream::out);

  size_t dim = nodes.get_numCols();

  outvtk << "# vtk DataFile Version 1.0\nresults.vtk  2D Unstructured Grid of Triangles\nASCII\n\nDATASET UNSTRUCTURED_GRID\nPOINTS "
           << nodes.get_numRows() << " float\n";
  for (size_t i = 0; i < nodes.get_numRows(); ++i) {
    for (size_t j = 0; j < nodes.get_numCols(); ++j) {
      outvtk << nodes(i, j) << " ";
    }
    outvtk << "0.0\n";
  }

  outvtk << "\nCELLS " << elements.get_numRows() << " " << elements.get_numRows() * 4 << "\n";
  for (size_t i = 0; i < elements.get_numRows(); ++i) {
    outvtk << dim + 1 << " ";
    for (size_t j = 0; j < elements.get_numCols(); ++j) {
      outvtk << elements(i, j) << " ";
    }
    outvtk << "\n";
  }

  outvtk << "\nCELL_TYPES " << elements.get_numRows() << "\n";
  for (size_t i = 0; i < elements.get_numRows(); ++i) {
    outvtk << 5 << "\n";
  }

  outvtk << "\nPOINT_DATA " << nodes.get_numRows() << "\n";
  outvtk << "\nVECTORS displacements float\n";
  for (size_t i = 0; i < displacements.get_numRows(); ++i) {
    for (size_t j = 0; j < displacements.get_numCols(); ++j) {
      outvtk << displacements(i, j) << " ";
    }
    outvtk << "0.0\n";
  }

  outvtk << "\nSCALARS summary float\nLOOKUP_TABLE default\n";
  for (size_t i = 0; i < displacements.get_numRows(); ++i) {
    float s = 0.f;
    for (size_t j = 0; j < displacements.get_numCols(); ++j) {
      s += displacements(i, j) * displacements(i, j);
    }
    outvtk << std::sqrt(s) << "\n";
  }
}

void MakeVTKfile3D(std::string output_vtk,
                   Matrix &nodes,
                   Matrix &elements,
                   Matrix &displacements) {
  CheckRunTime(__func__)
  std::fstream outvtk;
  outvtk.open(output_vtk, std::fstream::out);

  size_t dim = nodes.get_numCols();

  outvtk << "# vtk DataFile Version 1.0\nresults.vtk  3D Unstructured Grid of Triangles\nASCII\n\nDATASET UNSTRUCTURED_GRID\nPOINTS "
           << nodes.get_numRows() << " float\n";
  for (size_t i = 0; i < nodes.get_numRows(); ++i) {
    for (size_t j = 0; j < nodes.get_numCols(); ++j) {
      outvtk << nodes(i, j) << " ";
    }
    outvtk << "\n";
  }

  outvtk << "CELLS " << elements.get_numRows() << " " << elements.get_numRows() * 5 << "\n";
  for (size_t i = 0; i < elements.get_numRows(); ++i) {
    outvtk << dim + 1 << " ";
    for (size_t j = 0; j < elements.get_numCols(); ++j) {
      outvtk << elements(i, j) << " ";
    }
    outvtk << "\n";
  }

  outvtk << "CELL_TYPES " << elements.get_numRows() << "\n";
  for (size_t i = 0; i < elements.get_numRows(); ++i) {
    outvtk << 10 << "\n";
  }

  outvtk << "\nPOINT_DATA " << nodes.get_numRows() << "\n";

  outvtk << "\nVECTORS displacements float\n";
  for (size_t i = 0; i < displacements.get_numRows(); ++i) {
    for (size_t j = 0; j < displacements.get_numCols(); ++j) {
      outvtk << displacements(i, j) << " ";
    }
    outvtk << "\n";
  }

  outvtk << "\nSCALARS summary float\nLOOKUP_TABLE default\n";
  for (size_t i = 0; i < displacements.get_numRows(); ++i) {
    float s = 0.f;
    for (size_t j = 0; j < displacements.get_numCols(); ++j) {
      s += displacements(i, j) * displacements(i, j);
    }
    outvtk << std::sqrt(s) << "\n";
  }
}

void MakeVTKfile3D(std::string output_vtk,
                   std::vector<CPU_Matrix> &nodes,
                   std::vector<Element> &elements,
                   Matrix &displacements,
                   std::vector<CPU_Matrix> &Stress,
                   std::vector<float> &sigma_mises,
                   std::vector<CPU_Matrix> &Deformation,
                   std::vector<float> &epsilon_mises) {
  CheckRunTime(__func__)
  std::fstream outvtk;
  outvtk.open(output_vtk, std::fstream::out);
  outvtk << "# vtk DataFile Version 1.0\nresults.vtk  3D Unstructured Grid of Triangles\nASCII\n\nDATASET UNSTRUCTURED_GRID\nPOINTS "
           << nodes[0].get_numElements() << " float\n";
  for (size_t i = 0; i < nodes[0].get_numElements(); ++i) {
    outvtk << nodes[0][i] << " " << nodes[1][i] << " " << nodes[2][i] << "\n";
  }

  outvtk << "CELLS " << elements.size() << " " << elements.size() * 5 << "\n";
  for (size_t i = 0; i < elements.size(); ++i) {
    outvtk << 4 << " " << elements[i].nodesIds[0] << " " << elements[i].nodesIds[1] << " " << elements[i].nodesIds[2] << " " << elements[i].nodesIds[3] << "\n";
  }

  outvtk << "CELL_TYPES " << elements.size() << "\n";
  for (size_t i = 0; i < elements.size(); ++i) {
    outvtk << 10 << "\n";
  }

  outvtk << "\nPOINT_DATA " << nodes[0].get_numElements() << "\n";

  outvtk << "\nVECTORS displacements float\n";
  for (size_t i = 0; i < displacements.get_numElements() - 2; i += 3) {
    outvtk << displacements[i] << " " << displacements[i + 1] << " " << displacements[i + 2] << "\n";
  }

  outvtk << "\nSCALARS summary float\nLOOKUP_TABLE default\n";
  for (size_t i = 0; i < displacements.get_numElements() - 2; i += 3) {
    outvtk << std::sqrt(displacements[i] * displacements[i] + displacements[i + 1] * displacements[i + 1] + displacements[i + 2] * displacements[i + 2]) << "\n";
  }

  //    outvtk << "\nCELL_DATA " << elements.size() << "\n";
  //    outvtk << "VECTORS stress float\n";
  //    for (int i = 0; i < Stress.size(); i++) {
  //        outvtk << Stress[i][0] << " " << Stress[i][1] << " " << Stress[i][2] << " " << Stress[i][3] << " " << Stress[i][4] << " " << Stress[i][5] << "\n";
  //    }

//      outvtk << "\nSCALARS mises_stress float\nLOOKUP_TABLE default\n";
//      for (int i = 0; i < sigma_mises.size(); i++) {
//          outvtk << sigma_mises[i] << "\n";
//      }

  //    outvtk << "\nVECTORS deformation float\n";
  //    for (int i = 0; i < Deformation.size(); i++) {
  //        outvtk << Deformation[i][0] << " " << Deformation[i][1] << " " << Deformation[i][2] << " " << Deformation[i][3] << " " << Deformation[i][4] << " " << Deformation[i][5] << "\n";
  //    }

  //    outvtk << "\nSCALARS mises_deformation float\nLOOKUP_TABLE default\n";
  //    for (int i = 0; i < epsilon_mises.size(); i++) {
  //        outvtk << epsilon_mises[i] << "\n";
  //    }
}

