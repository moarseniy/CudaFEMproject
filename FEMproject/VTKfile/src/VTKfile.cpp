
#include "VTKfile.h"

void MakeVTKfile(std::string output_vtk,
                 MyArray nodesX,
                 MyArray nodesY,
                 MyArray nodesZ,
                 std::vector<Element> elements,
                 MyArray displacements,
                 std::vector<MyArray> Stress,
                 std::vector<float> sigma_mises,
                 std::vector<MyArray> Deformation,
                 std::vector<float> epsilon_mises) {
    #ifdef TOOLS_TIMER
        Timer timer(__func__);
    #endif
    fstream outvtk;
    outvtk.open(output_vtk, fstream::out);
    outvtk << "# vtk DataFile Version 1.0\nresults.vtk  3D Unstructured Grid of Triangles\nASCII\n\nDATASET UNSTRUCTURED_GRID\nPOINTS "
           << nodesX.get_size() << " float\n";
    for (int i = 0; i < nodesX.get_size(); i++) {
        outvtk << nodesX[i] << " " << nodesY[i] << " " << nodesZ[i] << "\n";
    }

    outvtk << "CELLS " << elements.size() << " " << elements.size() * 5 << "\n";
    for (int i = 0; i < elements.size(); i++) {
        outvtk << 4 << " " << elements[i].nodesIds[0] << " " << elements[i].nodesIds[1] << " " << elements[i].nodesIds[2] << " " << elements[i].nodesIds[3] << "\n";
    }

    outvtk << "CELL_TYPES " << elements.size() << "\n";
    for (int i = 0; i < elements.size(); i++) {
        outvtk << 10 << "\n";
    }

    outvtk << "\nPOINT_DATA " << nodesX.get_size() << "\n";

    outvtk << "\nVECTORS displacements float\n";
    for (int i = 0; i < displacements.get_size() - 2; i += 3) {
        outvtk << displacements[i] << " " << displacements[i + 1] << " " << displacements[i + 2] << "\n";
    }

    outvtk << "\nSCALARS summary float\nLOOKUP_TABLE default\n";
    for (int i = 0; i < displacements.get_size() - 2; i += 3) {
        outvtk << std::sqrt(displacements[i] * displacements[i] + displacements[i + 1] * displacements[i + 1] + displacements[i + 2] * displacements[i + 2]) << "\n";
    }

//    outvtk << "\nCELL_DATA " << elements.size() << "\n";
//    outvtk << "VECTORS stress float\n";
//    for (int i = 0; i < Stress.size(); i++) {
//        outvtk << Stress[i][0] << " " << Stress[i][1] << " " << Stress[i][2] << " " << Stress[i][3] << " " << Stress[i][4] << " " << Stress[i][5] << "\n";
//    }

//    outvtk << "\nSCALARS mises_stress float\nLOOKUP_TABLE default\n";
//    for (int i = 0; i < sigma_mises.size(); i++) {
//        outvtk << sigma_mises[i] << "\n";
//    }

//    outvtk << "\nVECTORS deformation float\n";
//    for (int i = 0; i < Deformation.size(); i++) {
//        outvtk << Deformation[i][0] << " " << Deformation[i][1] << " " << Deformation[i][2] << " " << Deformation[i][3] << " " << Deformation[i][4] << " " << Deformation[i][5] << "\n";
//    }

//    outvtk << "\nSCALARS mises_deformation float\nLOOKUP_TABLE default\n";
//    for (int i = 0; i < epsilon_mises.size(); i++) {
//        outvtk << epsilon_mises[i] << "\n";
//    }
}


void MakeVTKfileNew(std::string output_vtk,
                 MyArray nodesX,
                 MyArray nodesY,
                 MyArray nodesZ,
                 std::vector<ElementLight> elements,
                 MyArray displacements) {
    fstream outvtk;
    outvtk.open(output_vtk, fstream::out);
    outvtk << "# vtk DataFile Version 1.0\nresults.vtk  3D Unstructured Grid of Triangles\nASCII\n\nDATASET UNSTRUCTURED_GRID\nPOINTS "
           << nodesX.get_size() << " float\n";
    for (int i = 0; i < nodesX.get_size(); i++) {
        outvtk << nodesX[i] << " " << nodesY[i] << " " << nodesZ[i] << "\n";
    }

    outvtk << "CELLS " << elements.size() << " " << elements.size() * 5 << "\n";
    for (int i = 0; i < elements.size(); i++) {
        outvtk << 4 << " " << elements[i].nodesIds[0] << " " << elements[i].nodesIds[1] << " " << elements[i].nodesIds[2] << " " << elements[i].nodesIds[3] << "\n";
    }

    outvtk << "CELL_TYPES " << elements.size() << "\n";
    for (int i = 0; i < elements.size(); i++) {
        outvtk << 10 << "\n";
    }

    outvtk << "\nPOINT_DATA " << nodesX.get_size() << "\n";

    outvtk << "\nVECTORS displacements float\n";
    for (int i = 0; i < displacements.get_size() - 2; i += 3) {
        outvtk << displacements[i] << " " << displacements[i + 1] << " " << displacements[i + 2] << "\n";
    }

    outvtk << "\nSCALARS summary float\nLOOKUP_TABLE default\n";
    for (int i = 0; i < displacements.get_size() - 2; i += 3) {
        outvtk << std::sqrt(displacements[i] * displacements[i] + displacements[i + 1] * displacements[i + 1] + displacements[i + 2] * displacements[i + 2]) << "\n";
    }

//    outvtk << "\nCELL_DATA " << elements.size() << "\n";
//    outvtk << "VECTORS stress float\n";
//    for (int i = 0; i < Stress.size(); i++) {
//        outvtk << Stress[i][0] << " " << Stress[i][1] << " " << Stress[i][2] << " " << Stress[i][3] << " " << Stress[i][4] << " " << Stress[i][5] << "\n";
//    }

//    outvtk << "\nSCALARS mises_stress float\nLOOKUP_TABLE default\n";
//    for (int i = 0; i < sigma_mises.size(); i++) {
//        outvtk << sigma_mises[i] << "\n";
//    }

//    outvtk << "\nVECTORS deformation float\n";
//    for (int i = 0; i < Deformation.size(); i++) {
//        outvtk << Deformation[i][0] << " " << Deformation[i][1] << " " << Deformation[i][2] << " " << Deformation[i][3] << " " << Deformation[i][4] << " " << Deformation[i][5] << "\n";
//    }

//    outvtk << "\nSCALARS mises_deformation float\nLOOKUP_TABLE default\n";
//    for (int i = 0; i < epsilon_mises.size(); i++) {
//        outvtk << epsilon_mises[i] << "\n";
//    }
}

