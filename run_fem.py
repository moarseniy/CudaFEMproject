import os, json, sys, shutil

from mesh_parser_lib.fem_parser import run_parser
from postprocessing.postprocessing import run_postprocessing


if __name__ == '__main__':
    
    try:
        config_path = str(sys.argv[1])
        with open(config_path, 'r') as file_read:
            config_data = json.load(file_read)
            
        proj_dir = str(config_data['CudaFEMproject_path'])
        res_dir = str(config_data['results_path'])
        prep_meshes_dir = str(config_data['prepared_meshes_path'])
        task_list = config_data['task']
        dim = str(config_data['DIM'])
        youngModulus = str(config_data['youngModulus'])
        poissonRatio = str(config_data['poissonRatio'])

    except Exception as e:
        print('Something wrong with config file: ', e)
        sys.exit(0)

    for task_name in task_list:
        with_parser = True
        prep_mesh_task_dir = prep_meshes_dir + dim + 'D/' + task_name
        raw_mesh_dir = proj_dir + 'raw_meshes/' + dim + 'D/' + task_name + '.k'
        
        if not os.path.exists(raw_mesh_dir):
            print('raw mesh ' + task_name + '.k' + 'doesn\'t exist')
            sys.exit(0)
        if os.path.exists(prep_mesh_task_dir):
            with_parser = False
        else:
            os.mkdir(prep_mesh_task_dir)

        if not os.path.exists(res_dir + dim + 'D/' + task_name):
            os.mkdir(res_dir + dim + 'D/' + task_name)
        if os.path.exists(res_dir + dim + 'D/' + task_name + '/output/'):
            shutil.rmtree(res_dir + dim + 'D/' + task_name + '/output/')
            os.mkdir(res_dir + dim + 'D/' + task_name + '/output/')
        else:
            os.mkdir(res_dir + dim + 'D/' + task_name + '/output/')

        if with_parser:
            print('RUN parser for ' + task_name)
            run_parser(proj_dir, prep_mesh_task_dir, raw_mesh_dir, task_name, dim)
            print('FINISHED parser for ' + task_name)

        print('\nRUN calculation for ' + task_name)
        main_exe = 'build-FEMproject-Desktop_x86_windows_msvc2019_pe_64bit-Release/fem/main.exe'
        os.system(proj_dir + main_exe + ' ' + task_name + ' ' + proj_dir + ' ' + prep_meshes_dir + ' ' +
                  res_dir + ' ' + youngModulus + ' ' + poissonRatio)
        print('FINISHED calculation for ' + task_name)
        
        print('\nRUN Postprocessing for ' + task_name)
        run_postprocessing(res_dir + dim + 'D/' + task_name + '/output/')
        print('FINISHED Postprocessing for ' + task_name)
    

    
