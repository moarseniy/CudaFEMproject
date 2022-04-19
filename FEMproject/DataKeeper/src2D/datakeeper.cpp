
#include "datakeeper.h"

using namespace std;

// I always liked strtof since it lets you specify an end pointer. [ https://stackoverflow.com/a/57163016 ]
bool isFloat(const std::string& str)
{
    char* ptr;
    strtof(str.c_str(), &ptr);
    return (*ptr) == '\0';
}

void FEMdataKeeper::ParseFiles(float poissonRatio, float youngModulus) {
    CheckRunTime(__func__)

    fstream nodes_file, elements_file, loads_file, constraints_file, stress_file;

    nodes_file.open         (prep_mesh_dir + "/nodes.txt",         fstream::in);
    elements_file.open      (prep_mesh_dir + "/elements.txt",      fstream::in);
    loads_file.open         (prep_mesh_dir + "/loads.txt",         fstream::in);
    constraints_file.open   (prep_mesh_dir + "/constraints.txt",   fstream::in);
    stress_file.open        (prep_mesh_dir + "/stress.txt",        fstream::in);

    nodes_file          >> nodesCount;
    elements_file       >> elementsCount;
    constraints_file    >> constraintsCount;
    loads_file          >> loadsCount;
    stress_file         >> boundaryEdgesCount;

    AllocateDynamicMemory();

    D(0,0) = 1.0;			D(0, 1) = poissonRatio;	D(0, 2) = 0.0;
    D(1, 0) = poissonRatio;	D(1, 1) = 1.0; 			D(1, 2) = 0.0;
    D(2, 0) = 0.0;        	D(2, 1) = 0.0;        	D(2, 2) = (1.0 - poissonRatio) / 2.0;
    D.scale(youngModulus / (1.0 - pow(poissonRatio, 2.0)));

//    std::cout << "NODES" << std::endl;
    for (int i = 0; i < nodesCount; ++i) {
        nodes_file >> nodesX[i] >> nodesY[i];
//        std::cout << nodesX[i] << ' ' << nodesY[i] << std::endl;
    }

    for (int i = 0; i < elementsCount; ++i) {
        Element element;
        elements_file >> element.nodesIds[0] >> element.nodesIds[1] >> element.nodesIds[2];
        // Jacobian calucation
//        float X2 = nodesX[element.nodesIds[1]] - nodesX[element.nodesIds[0]];
//        float X3 = nodesX[element.nodesIds[2]] - nodesX[element.nodesIds[0]];
//        float Y2 = nodesY[element.nodesIds[1]] - nodesY[element.nodesIds[0]];
//        float Y3 = nodesY[element.nodesIds[2]] - nodesY[element.nodesIds[0]];
//        element.jacobian = 1.0 / (X3 * Y2 - X2 * Y3);
        elements.push_back(element);
    }

//    std::cout << "CONSTRAINTS" << std::endl;
    for (int i = 0; i < constraintsCount; ++i) {
        Constraint constraint;
        int type;
        constraints_file >> constraint.node >> type;
        constraint.type = static_cast<Constraint::Type>(type);
        constraints.push_back(constraint);
//        std::cout << constraint.node << ' ' << type << std::endl;
    }

//    std::cout << "LOADS" << std::endl;
    // ToDO: Add parsing of time-dependant functions like Ricker
    for (int i = 0; i < loadsCount; ++i) {
        int node; //float xampl, yampl;
        std::string xstr, ystr, str;

        //loads_file >> node >> xampl >> yampl;
        loads_file >> node;
        std::getline(loads_file, str);
        // https://stackoverflow.com/a/39080627
        std::cout << str << std::endl;
        std::regex wavelet_expression("(\\w+(\\((.*?)\\))|[+-]?(\\d+([.]\\d*)?([eE][+-]?\\d+)?|[.]\\d+([eE][+-]?\\d+)?))");
        smatch res;

        int loop_it = 0;
        while (std::regex_search(str, res, wavelet_expression)) {
            ++loop_it;
            if (loop_it == 1)      xstr = res[0];
            else if (loop_it == 2) ystr = res[0];
            str = res.suffix();
        }
        std::cout << "X:\n";
        if (isFloat(xstr)) {
            float xstr_value = std::stof(xstr);
            if (!(xstr_value == 0.0f)) {
                Load load;
                load.dof = node * 2 + 0;
                load.value = xstr_value;
                load.ampl = xstr_value;
                loads.push_back(load);
            }
        } else {             // if xstr is a time-dependent function, not constant value
            // parse string
            // https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c

            std::string wavelet_str = xstr.substr(0, xstr.find("("));
            xstr.erase(0, xstr.find("(") + 1);
            std::cout << "DEBUG : xstr = " << xstr << std::endl;
            std::cout << "DEBUG : wavelet_str = " << wavelet_str << std::endl;

            std::string freq_str = xstr.substr(0, xstr.find(","));
            xstr.erase(0, xstr.find(",") + 1);
            std::cout << "DEBUG : xstr = " << xstr << std::endl;
            std::cout << "DEBUG : freq_str = " << freq_str << std::endl;

            std::string timeshift_str = xstr.substr(0, xstr.find(","));
            xstr.erase(0, xstr.find(",") + 1);
            std::cout << "DEBUG : xstr = " << xstr << std::endl;
            std::cout << "DEBUG : timeshift_str = " << timeshift_str << std::endl;

            std::string ampl_str = xstr.substr(0, xstr.find(")"));
            //xstr.erase(0, xstr.find(")") + 1);
            std::cout << "DEBUG : xstr = " << xstr << std::endl;
            std::cout << "DEBUG : ampl_str = " << ampl_str << std::endl;

            Load load;
            load.dof   = node * 2 + 0;
            load.value = std::stof(ampl_str);
            load.ampl  = std::stof(ampl_str);

            // Try reading if you want more elegancy:
            // https://stackoverflow.com/questions/650162/why-the-switch-statement-cannot-be-applied-on-strings

            // https://stackoverflow.com/a/313990
            std::transform(wavelet_str.begin(), wavelet_str.end(), wavelet_str.begin(),
                                           [](unsigned char c){ return std::tolower(c); });
            if (wavelet_str == "ricker") {
                load.wavelet = &Load::Ricker;
            } // else if (wavelet_str == "berlage")

            std::cout << "DEBUG : wavelet_str = " << wavelet_str << std::endl;

            load.timeshift = std::stof(timeshift_str);
            load.freq      = std::stof(freq_str);
            std::cout << "load.value = " << load.value << "; load.timeshift = " << load.timeshift << "; load.freq = " << load.freq << std::endl;
            loads.push_back(load);

        }

        std::cout << "Y:\n";
        if (isFloat(ystr)) {
            float ystr_value = std::stof(ystr);
            if (!(ystr_value == 0.0f)) {
                Load load;
                load.dof = node * 2 + 1;
                load.value = ystr_value;
                load.ampl = ystr_value;
                loads.push_back(load);
            }
        } else {             // if ystr is a time-dependent function, not constant value
            // parse string
            // https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c

            std::string wavelet_str = ystr.substr(0, ystr.find("("));
            ystr.erase(0, ystr.find("(") + 1);
            std::cout << "DEBUG : wavelet_str = " << wavelet_str << std::endl;

            std::string freq_str = ystr.substr(0, ystr.find(","));
            ystr.erase(0, ystr.find(",") + 1);
            std::cout << "DEBUG : freq_str = " << freq_str << std::endl;

            std::string timeshift_str = ystr.substr(0, ystr.find(","));
            ystr.erase(0, ystr.find(",") + 1);
            std::cout << "DEBUG : timeshift_str = " << timeshift_str << std::endl;

            std::string ampl_str = ystr.substr(0, ystr.find(")"));
            //xstr.erase(0, ystr.find(")") + 1);
            std::cout << "DEBUG : ampl_str = " << ampl_str << std::endl;

            Load load;
            load.dof   = node * 2 + 1;
            load.value = std::stof(ampl_str);
            load.ampl  = std::stof(ampl_str);

            // Try reading if you want more elegancy:
            // https://stackoverflow.com/questions/650162/why-the-switch-statement-cannot-be-applied-on-strings

            // https://stackoverflow.com/a/313990
            std::transform(wavelet_str.begin(), wavelet_str.end(), wavelet_str.begin(),
                                           [](unsigned char c){ return std::tolower(c); });
            if (wavelet_str == "ricker") {
                load.wavelet = &Load::Ricker;
            } // else if (wavelet_str == "berlage")

            std::cout << "DEBUG : wavelet_str = " << wavelet_str << std::endl;

            load.timeshift = std::stof(timeshift_str);
            load.freq      = std::stof(freq_str);
            loads.push_back(load);

        }

//        if (!(yampl == 0.0)) {
//            Load load;
//            load.dof = node * 2 + 1;
//            // EDITED!!!! For Ricker
//            load.wavelet = &Load::Ricker;
//            load.value = yampl;
//            load.ampl = yampl;
//            load.freq = 15;
//            loads.push_back(load);
//        }

//        int node;
//        float x, y;
//        loads_file >> node >> x >> y;
//        loads[2 * node + 0] = x;
//        loads[2 * node + 1] = y;
//        std::cout << node << ' ' << x << ' ' << y << std::endl;
    }

    for (int i = 0; i < boundaryEdgesCount; ++i) {
        BoundaryEdge edge;
        float normal_x, normal_y;
        stress_file >> edge.node0 >> edge.node1 >> edge.adj_elem1 >> normal_x >> normal_y >> pressure[i];
        float normal_length = std::sqrt(normal_x * normal_x + normal_y * normal_y);
        edge.normal_x = normal_x / normal_length;
        edge.normal_y = normal_y / normal_length;
        boundary.push_back(edge);
    }
}
