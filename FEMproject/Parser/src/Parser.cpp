
#include "Parser.h"


void ParseFiles(std::string dir, std::string name) {
    fstream nodes_file, elements_file, loads_file, constraints_file;
    nodes_file.open(dir + name + "/nodes.txt", fstream::in);
    elements_file.open(dir + name + "/elements.txt", fstream::in);
    loads_file.open(dir + name + "/loads.txt", fstream::in);
    constraints_file.open(dir + name + "/constraints.txt", fstream::in);
}
