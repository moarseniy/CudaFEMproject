from mesh_parser_lib.parser_src.parser_lib import FileParser
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
    print('Write task dimension:')
    dim = 2#int(input())
    print('Write filename to parse mesh from:')
    filename = 'test_rect_pcg'#str(input())

    # UPDATE IT!!!
    raw_mesh_dir = 'C:/Users/mokin/Desktop/git/CudaFEMproject/raw_meshes/' + str(dim) + \
    'D/' + filename + '.k'
    prep_mesh_dir = 'C:/Users/mokin/Desktop/fem_stuff_test/prepared_meshes/' + str(dim) + \
    'D/' + filename
    proj_dir = 'C:/Users/mokin/Desktop/git/CudaFEMproject/'
    #
    
    run_parser(proj_dir, prep_mesh_dir, raw_mesh_dir, filename, dim)
