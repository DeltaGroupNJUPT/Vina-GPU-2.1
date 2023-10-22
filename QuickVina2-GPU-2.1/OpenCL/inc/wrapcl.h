#pragma once
#include "commonMacros.h"


void checkErr(cl_int err);
void SetupPlatform(cl_platform_id** platforms, cl_int* gpu_platform_id);
void SetupDevice(cl_platform_id* platforms, cl_device_id** devices, cl_int gpu_platform_id);
void SetupContext(cl_platform_id* platforms, cl_device_id* devices, cl_context* context, cl_uint num_device, cl_int gpu_platform_id);
void SetupQueue(cl_command_queue* queue, cl_context context, cl_device_id* devices);
void read_file(char** program_file, size_t* program_size, std::string file_path);
void read_n_file(char** program_file, size_t* program_size, std::string file_paths[], size_t num_file);
void SetupBuildProgramWithSource(cl_program program_cl, cl_program program_head, cl_device_id* devices, std::string include_path, std::string addtion);
void SaveProgramToBinary(cl_program program_cl, const char* file_name);
cl_program SetupBuildProgramWithBinary(cl_context context, cl_device_id* devices, const char* binary_file_name);
void SetupKernel(cl_kernel* kernels, cl_program* programs, size_t num_kernel, const char kernel_name[][50]);
void CreateDeviceBuffer(cl_mem* mem, cl_mem_flags flag, size_t size, cl_context context);
void SetKernelArg(cl_kernel kernel, cl_uint num, size_t size, const void* ptr);