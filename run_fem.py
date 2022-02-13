from mesh_parser_lib.fem_parser import run_parser
import os, json, sys

if __name__ == '__main__':
    try:
        config_path = str(sys.argv[1])
        with open(config_path, 'r') as file_read:
            config_data = json.load(file_read)
            
        proj_dir = str(config_data['CudaFEMproject_path'])
        res_dir = str(config_data['results_path'])
        prep_meshes_dir = str(config_data['prepared_meshes_path'])
        task = str(config_data['task'])
        dim = str(config_data['DIM'])
        parser = str(config_data['parser'])
        youngModulus = str(config_data['youngModulus'])
        poissonRatio = str(config_data['poissonRatio'])
        
        run_parser(proj_dir, prep_meshes_dir, task, dim)

        if not os.path.exists(res_dir + dim + 'D/' + task):
            os.mkdir(res_dir + dim + 'D/' + task)
        if not os.path.exists(res_dir + dim + 'D/' + task + '/output/'):
            os.mkdir(res_dir + dim + 'D/' + task + '/output/')

        main_exe = 'build-FEMproject-Desktop_x86_windows_msvc2019_pe_64bit-Release/fem/main.exe'
        os.system(proj_dir + main_exe + ' ' + task + ' ' + proj_dir + ' ' + prep_meshes_dir + ' ' +
                  res_dir + ' ' + youngModulus + ' ' + poissonRatio)
                
    except Exception as e:
        print('Something wrong with config file: ', e)
        sys.exit(0)

    
