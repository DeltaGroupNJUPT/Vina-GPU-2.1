#ifndef COMMON_MACROS_H
#define COMMON_MACROS_H
#define _CRT_SECURE_NO_WARNINGS
//#define CL_USE_DEPRECATED_OPENCL_3_0_APIS

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#ifdef __APPLE__
	#include "OpenCL/opencl.h"
#else
	#include "CL/opencl.h"
#endif

#ifdef WINDOWS
#include <Windows.h>
#else
#include <unistd.h>
#endif

#include "kernel2.h"
#define KERNEL_1_NUM_OF_KERNELS 1 
#define KERNEL_2_NUM_OF_KERNELS 2
#define NUM_OF_FILES 5
#define NUM_OF_FILES_KERNEL_1 2
#define NUM_OF_FILES_KERNEL_2 5
#endif

