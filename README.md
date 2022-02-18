# CudaFEMproject
Finite Element Method (3D/2D) + CUDA + CMake implementation
* Arseniy Mokin moarseniy@yandex.ru
* Grigory Sabinin gvsabinin@gmail.com

## How to use?
0. [Description](#Description)
1. [Required software](#Required-software)
2. [Run project](#Run-project)
3. [Postprocessing](#Postprocessing)

![vmY5Emnim5I](https://user-images.githubusercontent.com/44135971/144502755-0cba54a7-e0bf-4187-9dc5-5f43491439fd.jpg)
![image](https://user-images.githubusercontent.com/44135971/151878668-b4dedb38-d3af-46b8-bf9e-9f80ae676071.png)
![9BjXtg9EDbE](https://user-images.githubusercontent.com/44135971/144608765-c2a25a08-a927-4e12-89b3-3eed081d221b.jpg)

## Description
Project

## Required software
* gcc (C/C++)
* CMake
* Python 3.7+
* QtCreator (optional)
* CUDA Toolkit (optinal)

## Run project
* git clone https://github.com/moarseniy/CudaFEMproject.git
* Open CMakeLists.txt file contained in repository through Qt Creator
* Run project in Release mode
* Open and change run_config in configs directory

'''javascript
{
	"CudaFEMproject_path": "C:/Users/mokin/Desktop/git/CudaFEMproject/",
	"results_path": "C:/Users/mokin/Desktop/fem_stuff/fem_results/",
	"prepared_meshes_path": "C:/Users/mokin/Desktop/fem_stuff/prepared_meshes/",
	"task": ["9task_3", "9task_4", "9task_5"],
	"DIM": 2,
	"parser": true,
	"youngModulus": 0.25,
	"poissonRatio": 2e+07
}
'''

* Open terminal in CudaFEMproject and start **python run_fem.py C:/Users/mokin/Desktop/git/CudaFEMproject/configs/run_config.json**

## Postprocessing

coming soon




