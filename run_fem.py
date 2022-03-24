import os, json, sys, shutil

from mesh_parser_lib.fem_parser import run_parser
from postprocessing.postprocessing import run_postprocessing


if __name__ == '__main__':
    
    try:
        config_path = str(sys.argv[1])
        with open(config_path, 'r') as file_read:
            config_data = json.load(file_read)
            
        proj_dir = str(config_data['CudaFEMproject_path'])
        work_dir = str(config_data['working_path'])
        task_list = config_data['task']
        dim = str(config_data['DIM'])
        youngModulus = str(config_data['youngModulus'])
        poissonRatio = str(config_data['poissonRatio'])

    except Exception as e:
        print('Something wrong with config file: ', e)
        sys.exit(0)

    for task_name in task_list:
        with_parser = True
        
        raw_mesh_dir = proj_dir + 'raw_meshes/' + dim + 'D/' + task_name + '.k'
        prep_mesh_task_dir = work_dir + 'prepared_meshes/' + dim + 'D/' + task_name
        results_dir = work_dir + 'fem_results/' + dim + 'D/' + task_name

        #raw_meshes
        if not os.path.exists(raw_mesh_dir):
            print('raw mesh ' + task_name + '.k' + 'doesn\'t exist')
            sys.exit(0)

        #prep_meshes
        if os.path.exists(prep_mesh_task_dir):
            with_parser = False
        else:
            if not os.path.exists(work_dir + 'prepared_meshes/'):
                os.mkdir(work_dir + 'prepared_meshes/')
            if not os.path.exists(work_dir + 'prepared_meshes/' + dim + 'D/'):
                os.mkdir(work_dir + 'prepared_meshes/' + dim + 'D/')
            os.mkdir(prep_mesh_task_dir)

        #final_results
        if os.path.exists(results_dir + '/output/'):
            shutil.rmtree(results_dir + '/output/')
            os.mkdir(results_dir + '/output/')
        else:
            if not os.path.exists(work_dir + 'fem_results/'):
                os.mkdir(work_dir + 'fem_results/')
            if not os.path.exists(work_dir + 'fem_results/' + dim + 'D/'):
                os.mkdir(work_dir + 'fem_results/' + dim + 'D/')
            if not os.path.exists(results_dir):
                os.mkdir(results_dir)
            os.mkdir(results_dir + '/output/')

        ########################################
            
        if with_parser:
            print('RUN parser for ' + task_name)
            run_parser(proj_dir, prep_mesh_task_dir, raw_mesh_dir, task_name, dim)
            print('FINISHED parser for ' + task_name)

        print('\nRUN calculation for ' + task_name)
        main_exe = 'build-FEMproject-Desktop_x86_windows_msvc2019_pe_64bit-Release/fem/main.exe'
        os.system(proj_dir + main_exe + ' ' + task_name + ' ' + proj_dir + ' ' + prep_mesh_task_dir + ' ' +
                  results_dir + ' ' + poissonRatio + ' ' + youngModulus)
        print('FINISHED calculation for ' + task_name)
        
        print('\nRUN Postprocessing for ' + task_name)
        run_postprocessing(results_dir + '/output/')
        print('FINISHED Postprocessing for ' + task_name)
    

    
