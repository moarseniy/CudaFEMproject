from mesh_parser_lib.parser import FileParser
import os

def run_parser(proj_dir, prepared_mesh_dir, raw_mesh_dir, filename, dim):

    filename += '.k'
    task_name = os.path.splitext(filename)[0]

    parser = FileParser(filename, prepared_mesh_dir, raw_mesh_dir, task_name, int(dim))
    #parser.make_directory()

    parser.parse_nodes()
    parser.prepare_nodes()

    parser.parse_loads()
    
    parser.parse_elements()
    parser.prepare_elements()

    parser.parse_constraints_and_sets()
    parser.prepare_constraints()

    parser.prepare_loads()


if __name__ == '__main__':
    # UPDATE IT!!!
    print('Write task dimension:')
    dim = int(input())
    print('Write filename to parse mesh from:')
    filename = str(input())                   
    run_parser(proj_dir, prep_mesh_dir, filename, dim)
