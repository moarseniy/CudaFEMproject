from mesh_parser_lib.file_parser import FileParser
import os

def run_parser(filename, dim):
    
    raw_mesh_dir = 'raw_meshes/' + str(dim) + 'D'   # Folder with raw meshes
    prepared_mesh_dir = 'prepared_meshes/' + str(dim) + 'D'  # Folder for prepared meshes

    if not os.path.splitext(filename)[1]:
        filename += '.k'
    task_name = os.path.splitext(filename)[0]         # Folder for results

    parser = FileParser(filename, prepared_mesh_dir, raw_mesh_dir, task_name, dim)
    parser.make_directory()

    parser.parse_nodes()
    parser.prepare_nodes()

    parser.parse_loads()
    
    parser.parse_elements()
    parser.prepare_elements()

    parser.parse_constraints_and_sets()
    parser.prepare_constraints()

    parser.prepare_loads()


if __name__ == '__main__':
    print('Write task dimension:')
    dim = int(input())
    print('Write filename to parse mesh from:')
    filename = str(input())                     
    run_parser(filename, dim)
