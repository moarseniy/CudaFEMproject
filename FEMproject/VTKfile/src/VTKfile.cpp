
#include <VTKfile.h>


void writeSEGY(std::string path, Matrix &v) {
  CSegyOut out(path);
  out.open();

//  out.setSegyFormat();
//  out.setSegyID();
//  out.setSegyTrace();
//  out.setTraceNumber();
//  out.setSampleInterval();
//  out.setSampleNumber();

  out.write();
  out.close();
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

