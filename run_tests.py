import os, json, sys, shutil
import numpy as np

def run_tests(proj_dir, work_dir, task_list):
    for task_name in task_list:
        try:
            config_path = proj_dir + 'testing/' + task_name + '_test/run_config.json'
            with open(config_path, 'r') as file_read:
                config_data = json.load(file_read)

            dim = str(config_data['DIM'])
            youngModulus = str(config_data['youngModulus'])
            poissonRatio = str(config_data['poissonRatio'])

        except Exception as e:
            print('Something wrong with config file: ', e)
            sys.exit(0)
        
        raw_mesh_dir = proj_dir + 'raw_meshes/' + dim + 'D/' + task_name + '.k'
        prep_mesh_task_dir = proj_dir + 'testing/' + task_name + '_test/' + task_name + '_mesh/'
        results_dir = work_dir + 'fem_results/' + dim + 'D/' + task_name

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
            
        main_exe = 'build-FEMproject-Desktop_x86_windows_msvc2019_pe_64bit-Release/testing/testing.exe'
        os.system(proj_dir + main_exe + ' ' + task_name + ' ' + proj_dir + ' ' + prep_mesh_task_dir + ' ' +
                  results_dir + ' ' + poissonRatio + ' ' + youngModulus)

        eps = 1e-30
        x_res = []
        y_res = []
        z_res = []
        with open(results_dir + '/results.vtk', 'r') as res_file:
            for line in res_file:
                if 'displacements' in line:
                    line = res_file.readline()
                    while line != '\n':
                        x_res.append(float(line.split(' ')[0]))
                        y_res.append(float(line.split(' ')[1]))
                        z_res.append(float(line.split(' ')[2]))
                        line = res_file.readline()
        x_test = []
        y_test = []
        z_test = []
        with open(proj_dir + 'testing/' + task_name + '_test/results.vtk', 'r') as res_file:
            for line in res_file:
                if 'displacements' in line:
                    line = res_file.readline()
                    while line != '\n':
                        x_test.append(float(line.split(' ')[0]))
                        y_test.append(float(line.split(' ')[1]))
                        z_test.append(float(line.split(' ')[2]))
                        line = res_file.readline()
        if np.linalg.norm(np.array(x_res) - np.array(x_test)) < eps and \
            np.linalg.norm(np.array(y_res) - np.array(y_test)) < eps and \
            np.linalg.norm(np.array(z_res) - np.array(z_test)) < eps:
                print('Test for ' + task_name + ' completed successfuly')

