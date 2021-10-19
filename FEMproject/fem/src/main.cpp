
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <time.h>
#include <direct.h>     // for _mkdir (Visual Studio specific)

#include "Tools.h"
#include "Linal2.h"
//#include <cmath>
#include "femfunc.h"
#include "VTKfile.h"
#include "init.h"
#include "datakeeper.h"

using namespace std;

void test(void);

int main(void) {
    CheckRunTime(__func__)

    std::string name = "bulk_small";

    std::string project_directory = "C:/Users/mexika/Documents/Qt_code/CudaFEMproject/";
    std::string mesh_directory = project_directory + "prepared_meshes/";
    std::string results_directory = project_directory + "final_results/" + name + "/";
    _mkdir(results_directory.c_str());
    std::string output_vtk = results_directory + "results.vtk";
    std::string output_results = results_directory + "output.txt";

    float poissonRatio = 0.3, youngModulus = 2e+11;

    FEMdataKeeper FEMdata;
    FEMdata.ParseFiles(mesh_directory, name, poissonRatio, youngModulus);
    FEMdata.ShowInfo();

    CalculateFiniteElementMethod(FEMdata);

    MakeResults(FEMdata, output_vtk);

    return 0;
}

void test(void)
{
    SparseMatrixCOO mat(5);
    mat.add_value(2,1,0);
    mat.add_value(1,4,0);
    mat.add_value(4,4,0);
    mat.add_value(1,1,0);
    mat.add_value(3,4,0);

    mat.Show();

    SparseMatrixCOO mat2(2);
    //mat2 = mat.DeleteZeros();
    //mat2.Show();
    mat2 = mat.DeleteZeros();
    std::cout << "mat2.ptr = " << mat2.get_ptr() << std::endl;
    mat2.Show();
    mat.DeleteZerosOrder();
    mat.Show();
}
