# CMakeLists.txt for building Vina-GPU 2.1 on Windows with VS2022 + Boost 1.83 + CUDA 13.2 OpenCL
# Custom build script (not from upstream — upstream ships only a Linux Makefile)
# Author: written for Dhruvik N. Chauhan's kidney_lab project, Day 9 GPU investigation
# Target: AutoDock-Vina-GPU-2.1 binary on Blackwell sm_120 (RTX 5050)

cmake_minimum_required(VERSION 3.20)

project(VinaGPU LANGUAGES CXX)

# C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default Release build
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# ---------------------------------------------------------------------------
# Boost location hint — installer prebuilt at C:/local/boost_1_83_0
# CMake 4.x removed FindBoost.cmake; we must use BoostConfig.cmake (config mode)
# ---------------------------------------------------------------------------
set(BOOST_ROOT "C:/local/boost_1_83_0" CACHE PATH "Boost root directory")
set(Boost_ROOT "C:/local/boost_1_83_0" CACHE PATH "Boost root (CMake 4.x case)")
set(Boost_DIR "C:/local/boost_1_83_0/lib64-msvc-14.3/cmake/Boost-1.83.0" CACHE PATH "BoostConfig.cmake location")
list(APPEND CMAKE_PREFIX_PATH "C:/local/boost_1_83_0/lib64-msvc-14.3/cmake")

set(Boost_USE_STATIC_LIBS OFF)         # use the DLL versions
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

# Set CMake policies for Boost handling on CMake 4.x
if(POLICY CMP0144)
    cmake_policy(SET CMP0144 NEW)      # honor <PackageName>_ROOT variables
endif()
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)      # use BoostConfig.cmake (config mode) only
endif()

# Force config mode (CONFIG keyword); avoids legacy FindBoost lookup
# Include serialization — pulled in transitively by boost/thread internals
find_package(Boost 1.77 REQUIRED CONFIG COMPONENTS program_options system filesystem thread serialization)
message(STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost libraries: ${Boost_LIBRARIES}")
message(STATUS "Boost version: ${Boost_VERSION}")

# ---------------------------------------------------------------------------
# OpenCL — comes from CUDA 13.2 installation on Windows
# ---------------------------------------------------------------------------
set(CUDA_ROOT "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v13.2" CACHE PATH "CUDA root for OpenCL")
set(OpenCL_INCLUDE_DIR "${CUDA_ROOT}/include" CACHE PATH "OpenCL include directory")
set(OpenCL_LIBRARY "${CUDA_ROOT}/lib/x64/OpenCL.lib" CACHE FILEPATH "OpenCL library")

if(NOT EXISTS "${OpenCL_LIBRARY}")
    message(FATAL_ERROR "OpenCL.lib not found at ${OpenCL_LIBRARY}")
endif()
if(NOT EXISTS "${OpenCL_INCLUDE_DIR}/CL/cl.h")
    message(FATAL_ERROR "cl.h not found at ${OpenCL_INCLUDE_DIR}/CL/cl.h")
endif()

message(STATUS "OpenCL include: ${OpenCL_INCLUDE_DIR}")
message(STATUS "OpenCL library: ${OpenCL_LIBRARY}")

# ---------------------------------------------------------------------------
# Vina-GPU 2.1 source files
# ---------------------------------------------------------------------------
set(VINA_ROOT "${CMAKE_SOURCE_DIR}/AutoDock-Vina-GPU-2.1")

# Glob all .cpp files in lib/ (75 files)
file(GLOB VINA_LIB_SRC "${VINA_ROOT}/lib/*.cpp")
list(LENGTH VINA_LIB_SRC VINA_LIB_COUNT)
message(STATUS "Found ${VINA_LIB_COUNT} source files in lib/")

# wrapcl.cpp and main.cpp (per the upstream Makefile)
set(WRAPCL_SRC "${VINA_ROOT}/OpenCL/src/wrapcl.cpp")
set(MAIN_SRC "${VINA_ROOT}/main/main.cpp")

if(NOT EXISTS "${WRAPCL_SRC}")
    message(FATAL_ERROR "wrapcl.cpp not found at ${WRAPCL_SRC}")
endif()
if(NOT EXISTS "${MAIN_SRC}")
    message(FATAL_ERROR "main.cpp not found at ${MAIN_SRC}")
endif()

# ---------------------------------------------------------------------------
# Executable target
# ---------------------------------------------------------------------------
add_executable(AutoDock-Vina-GPU-2-1
    ${MAIN_SRC}
    ${WRAPCL_SRC}
    ${VINA_LIB_SRC}
)

target_include_directories(AutoDock-Vina-GPU-2-1 PRIVATE
    "${VINA_ROOT}/lib"
    "${VINA_ROOT}/OpenCL/inc"
    "${Boost_INCLUDE_DIRS}"
    "${OpenCL_INCLUDE_DIR}"
)

# Preprocessor definitions (matching upstream Makefile + README Windows recommendations)
target_compile_definitions(AutoDock-Vina-GPU-2-1 PRIVATE
    NVIDIA_PLATFORM                       # tells Vina-GPU we're on NVIDIA OpenCL
    OPENCL_2_0                            # OpenCL 3.0 (README recommends over 2.0)
    WINDOWS                               # Windows-specific code paths
    SMALL_BOX                             # default docking box size (can switch to LARGE_BOX later)
    BOOST_TIMER_ENABLE_DEPRECATED         # Boost 1.77+ requires this for the deprecated timer API
    BUILD_KERNEL_FROM_SOURCE              # compile .cl kernels at runtime from source (vs prebuilt .bin)
    NDEBUG                                # release-mode asserts off
    _SCL_SECURE_NO_WARNINGS               # quiet MSVC's iterator warnings
    _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS                # Boost uses some deprecated C++17 features
    _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING   # Vina-GPU main.cpp uses <experimental/filesystem>
    BOOST_ALLOW_DEPRECATED_HEADERS        # tolerate older boost headers
    BOOST_ALL_DYN_LINK                    # force Boost auto-link to use DLLs (matches Boost_USE_STATIC_LIBS OFF)
    CL_TARGET_OPENCL_VERSION=300          # explicitly target OpenCL 3.0 (silences cl_version.h info msg)
    # NOTE: _CRT_SECURE_NO_WARNINGS is defined inside Vina-GPU's commonMacros.h, so we don't redefine
)

target_link_libraries(AutoDock-Vina-GPU-2-1 PRIVATE
    Boost::program_options
    Boost::system
    Boost::filesystem
    Boost::thread
    Boost::serialization
    "${OpenCL_LIBRARY}"
)

# ---------------------------------------------------------------------------
# Windows / MSVC specific flags
# ---------------------------------------------------------------------------
if(MSVC)
    # /Zc:preprocessor — standards-conforming preprocessor (required by CUDA 13.x CCCL headers)
    # /bigobj — allow large object files (Vina-GPU has templates that generate huge .obj)
    # /MP — parallel compilation
    # /EHsc — standard C++ exception handling
    target_compile_options(AutoDock-Vina-GPU-2-1 PRIVATE
        /Zc:preprocessor
        /bigobj
        /MP
        /EHsc
        /wd4267    # silence size_t->int conversion warnings (Vina-GPU has many)
        /wd4244    # silence float->double conversion warnings
        /wd4996    # silence deprecated function warnings
    )
endif()

# Set output directory to project root for easier testing
set_target_properties(AutoDock-Vina-GPU-2-1 PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/build/bin"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_SOURCE_DIR}/build/bin"
)

message(STATUS "")
message(STATUS "============================================================")
message(STATUS "Vina-GPU 2.1 build configuration complete")
message(STATUS "  Boost     : ${Boost_VERSION} (${BOOST_ROOT})")
message(STATUS "  OpenCL    : ${OpenCL_LIBRARY}")
message(STATUS "  Source    : ${VINA_LIB_COUNT} files in lib/ + wrapcl.cpp + main.cpp")
message(STATUS "  Output    : build/bin/AutoDock-Vina-GPU-2-1.exe")
message(STATUS "============================================================")
