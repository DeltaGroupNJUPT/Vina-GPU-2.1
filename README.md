# Vina-GPU 2.1
Vina-GPU 2.1 further improves the virtual screening runtime and accuracy with the noval RILC-BFGS and GCS mthods based on Vina-GPU 2.0. 
Vina-GPU 2.1 includes AutoDock-Vina-GPU 2.1, QuickVina 2-GPU 2.1 and QuickVina-W-GPU 2.1.
![Vina-GPU2 1-arch](https://github.com/DeltaGroupNJUPT/Vina-GPU-2.1/assets/48940269/3b42ed59-01ce-449a-b203-deea1f0d0a36)



## Virtual Screening Results

* Runtime comparison of Vina-GPU 2.1 on Drugbank library (partial)
![runtime_all_drugbank](https://github.com/DeltaGroupNJUPT/Vina-GPU-2.1/assets/48940269/d728fee5-d4a6-4a16-bbde-cec06a81e38d)

* Accuracy comparison of Vina-GPU 2.1 on Drugbank library (partial)
![vs_accuracy_all_drugbank](https://github.com/DeltaGroupNJUPT/Vina-GPU-2.1/assets/48940269/eacffd1d-cb2a-40d9-9a74-4e8d071aac7b)

## Compiling and Running
### Windows
#### Build from source file
>Visual Studio 2019 is recommended for build Vina-GPU 2.1 from source
1. install [boost library](https://www.boost.org/) (current version is 1.77.0)
2. install [CUDA Toolkit](https://developer.nvidia.com/cuda-downloads) (current version: v12.2) if you are using NVIDIA GPU cards

    **Note**: the OpenCL library can be found in CUDA installation path for NVIDIA or in the driver installation path for AMD

3. add `$(ONE_OF_VINA_GPU_2_1_METHODS)/lib` `$(ONE_OF_VINA_GPU_2_1_METHODS)/OpenCL/inc` `$(YOUR_BOOST_LIBRARY_PATH)` `$(YOUR_BOOST_LIBRARY_PATH)/boost` `$(YOUR_CUDA_TOOLKIT_LIBRARY_PATH)/CUDA/v12.2/include` in the include directories
4. add `$(YOUR_BOOST_LIBRARY_PATH)/stage/lib` `$(YOUR_CUDA_TOOLKIT_PATH)/CUDA/lib/x64`in the addtional library 
5. add `OpenCL.lib` in the additional dependencies 
6. add `--config=$(ONE_OF_VINA_GPU_2_1_METHODS)/input_file_example/2bm2_config.txt` in the command arguments
7.  add `NVIDIA_PLATFORM` `OPENCL_3_0` `WINDOWS` in the preprocessor definitions if necessary
8. if you want to compile the binary kernel file on the fly, add `BUILD_KERNEL_FROM_SOURCE` in the preprocessor definitions
9. build & run
**Note**: ensure the line ending are CLRF

### Linux
**Note**: At least 8M stack size is needed. To change the stack size, use `ulimit -s 8192`.
1. install [boost library](https://www.boost.org/) (current version is 1.77.0)
2. install [CUDA Toolkit](https://developer.nvidia.com/cuda-downloads) (current version: v12.2) if you are using NVIDIA GPU cards

    **Note**: OpenCL library can be usually in `/usr/local/cuda` (for NVIDIA GPU cards)
3. `cd` into one of the three methods of Vina-GPU 2.1 (`$(ONE_OF_VINA_GPU_2_1_METHODS)`)
4. change the `BOOST_LIB_PATH` and `OPENCL_LIB_PATH` accordingly in `Makefile`
5. set GPU platform `GPU_PLATFORM` and OpenCL version `OPENCL_VERSION` in `Makefile`. some options are given below:

    **Note**: `-DOPENCL_3_0` is highly recommended in Linux, please avoid using `-OPENCL_1_2` in the Makefile setting. To check the OpenCL version on a given platform, use `clinfo`.
    |Macros|Options|Descriptions|
    |--|--|--|	
    |GPU_PLATFORM|-DNVIDIA_PLATFORM / -DAMD_PLATFORM|NVIDIA / AMD GPU platform
    |  OPENCL_VERSION | -DOPENCL_3_0 / -OPENCL_2_0|OpenCL version 2.1 / 2.0
    
6. type `make clean` and `make source` to build `$(ONE_OF_VINA_GPU_2_1_METHODS)` that compile the kernel files on the fly (this would take some time at the first use)
7. after a successful compiling, `$(ONE_OF_VINA_GPU_2_1_METHODS)` can be seen in the directory 
8. change `--opencl_binary_path` in the `./input_file_example/2bm2_config.txt` accordingly and type `$(ONE_OF_VINA_GPU_2_1_METHODS) --config ./input_file_example/2bm2_config.txt` to run one of the Vina-GPU 2.1 method
9. once you successfully run `$(ONE_OF_VINA_GPU_2_1_METHODS)`, its runtime can be further reduced by typing `make clean` and `make` to build it without compiling kernel files (but make sure the `Kernel1_Opt.bin` file and `Kernel2_Opt.bin` file is located in the dir specified by `--opencl_binary_path`)
10. other compile options: 

|Options| Description|
|--|--|
| -g | debug|
|-DTIME_ANALYSIS|output runtime analysis in `gpu_runtime.log`|
|-DDISPLAY_ADDITION_INFO|print addition information

## Enlarge the docking box
The docking box now can be enlarged by

1. Change `-DSMALL_BOX` into `-DLARGE_BOX` in `Makefile`
2. Type `make source` and `$(ONE_OF_VINA_GPU_2_1_METHODS) --config ./input_file_example/2bm2_config.txt`
3. Once the tutorial docking is finished, type `make clearn` and `make`
4. Now you can enlarge the docking box `--size_x/y/z` accordingly (see **Limitation** below)

## Structure Optimization
|Optimization|Methods|Reference|
|--|--|--|	
|Receptor Preparation|cross-docking| [origin paper](https://doi.org/10.1016/j.bmc.2022.116686)
|Binding Pocket Prediction|COACH-D|[origin paper](https://doi.org/10.1093/nar/gky439)
|Ligand Optimization| Gypsum-DL | [origin paper](https://doi.org/10.1186/s13321-019-0358-3)

## Usage
|Arguments| Description|Default value
|--|--|--|
|--config | the config file (in .txt format) that contains all the following arguments for the convenience of use| no default
| --receptor | the recrptor file (in .pdbqt format)| no default
|--ligand_directory| this path specifies the directory of all the input ligands(in .pdbqt format) | no default
|--output_directory| this path specifies the directory of the output ligands | no default
|--lbfgs| `--rilc_bfgs 0` turns off the RILC-BFGS, `--rilc_bfgs 1` turns on the RILC-BFGS` | `--rilc_bfgs 1`
|--thread| the scale of parallelism | 5000 for quickvina2-gpu 2.1, 8000 for others
|--search_depth| the number of searching iterations in each docking lane| heuristically determined
|--center_x/y/z|the center of searching box in the receptor|no default
|--size_x/y/z|the volume of the searching box|no default 
|--opencl_binary_path| this path specifies the directory of the kernel path | no default

## Limitation
|Arguments| Description|Limitation
|--|--|--|
|--thread| the scale of parallelism (docking lanes)| preferably less than 10000
|--size_x/y/z|the volume of the searching box |less than 100/100/100 for AutoDock-Vina-GPU 2.1 and 70/70/70 for other two variants

## Citation
* Ding J, Tang S, Mei Z, et al. Vina-GPU 2.0: Further Accelerating AutoDock Vina and Its Derivatives with Graphics Processing Units[J]. Journal of Chemical Information and Modeling, 2023, 63(7): 1982-1998.
* Tang, Shidi et al. “Accelerating AutoDock Vina with GPUs.” Molecules (Basel, Switzerland) vol. 27,9 3041. 9 May. 2022, doi:10.3390/molecules27093041
* Trott, Oleg, and Arthur J. Olson. "AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading." Journal of computational chemistry 31.2 (2010): 455-461.
* Hassan, N. M. , et al. "Protein-Ligand Blind Docking Using QuickVina-W With Inter-Process Spatio-Temporal Integration." Scientific Reports 7.1(2017):15451.
* Amr Alhossary, Stephanus Daniel Handoko, Yuguang Mu, and Chee-Keong Kwoh. "Fast, accurate, and reliable molecular docking with QuickVina 2. " Bioinformatics (2015): 2214–2216.
