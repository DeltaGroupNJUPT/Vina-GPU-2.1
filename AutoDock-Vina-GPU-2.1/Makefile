# Need to be modified according to different users
WORK_DIR=/home/shidi/Vina-GPU-2.1/AutoDock-Vina-GPU-2.1
BOOST_LIB_PATH=../../boost_1_84_0
OPENCL_LIB_PATH=/usr/local/cuda
OPENCL_VERSION=-DOPENCL_2_0
GPU_PLATFORM=-DNVIDIA_PLATFORM
DOCKING_BOX_SIZE=-DSMALL_BOX

# Should not be modified
BOOST_INC_PATH=-I$(BOOST_LIB_PATH) -I$(BOOST_LIB_PATH)/boost 
VINA_GPU_INC_PATH=-I$(WORK_DIR)/lib -I$(WORK_DIR)/OpenCL/inc 
OPENCL_INC_PATH=-I$(OPENCL_LIB_PATH)/include
LIB1=-lboost_program_options -lboost_system -lboost_filesystem -lOpenCL
LIB2=-lstdc++ -lstdc++fs
LIB3=-lm -lpthread
LIB_PATH=-L$(BOOST_LIB_PATH)/stage/lib -L$(OPENCL_LIB_PATH)/lib64
SRC=./lib/*.cpp ./OpenCL/src/wrapcl.cpp $(BOOST_LIB_PATH)/libs/thread/src/pthread/thread.cpp $(BOOST_LIB_PATH)/libs/thread/src/pthread/once.cpp #../boost_1_77_0/boost/filesystem/path.hpMACRO=-DAMD_PLATFORM -DDISPLAY_SUCCESS -DDISPLAY_ADDITION_INFO
MACRO=$(OPENCL_VERSION) $(GPU_PLATFORM) $(DOCKING_BOX_SIZE) -DBOOST_TIMER_ENABLE_DEPRECATED #-DDISPLAY_SUCCESS -DDISPLAY_ADDITION_INFO
all:out
out:./main/main.cpp
	gcc -o AutoDock-Vina-GPU-2-1 $(BOOST_INC_PATH) $(VINA_GPU_INC_PATH) $(OPENCL_INC_PATH) ./main/main.cpp -O3 $(SRC) $(LIB1) $(LIB2) $(LIB3) $(LIB_PATH) $(MACRO) $(OPTION) -DNDEBUG 
source:./main/main.cpp
	gcc -o AutoDock-Vina-GPU-2-1 $(BOOST_INC_PATH) $(VINA_GPU_INC_PATH) $(OPENCL_INC_PATH) ./main/main.cpp -O3 $(SRC) $(LIB1) $(LIB2) $(LIB3) $(LIB_PATH) $(MACRO) $(OPTION) -DNDEBUG -DBUILD_KERNEL_FROM_SOURCE 
debug:./main/main.cpp
	gcc -o AutoDock-Vina-GPU-2-1 $(BOOST_INC_PATH) $(VINA_GPU_INC_PATH) $(OPENCL_INC_PATH) ./main/main.cpp -g $(SRC) $(LIB1) $(LIB2) $(LIB3) $(LIB_PATH) $(MACRO) $(OPTION) -DBUILD_KERNEL_FROM_SOURCE 
clean:
	rm AutoDock-Vina-GPU-2-1
