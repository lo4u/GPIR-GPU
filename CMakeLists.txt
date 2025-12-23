cmake_minimum_required(VERSION 3.12)
project(spiral VERSION 0.1 LANGUAGES CXX CUDA)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_ARCHITECTURES 86)

add_compile_options(-Wno-deprecated-declarations)
set(CMAKE_CXX_FLAGS "-g -O3 -MD -march=native")

option(NOAVX512 "Disable AVX-512" OFF)
option(NOAVX2 "Disable AVX2" OFF)
if (NOAVX512)
message("Disabling AVX-512")
add_compile_options(-mno-avx512f)
endif()
if (NOAVX2)
message("Disabling AVX2")
add_compile_options(-mno-avx2)
endif()
if (NOCRT)
message("Disabling CRT")
add_compile_options(-DNOCRT=1)
endif()

# main target
# 设置可执行文件的输出目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
set(SRC_FILES src/spiral.cpp src/core.cpp src/constants.cpp src/poly.cpp src/util.cpp src/client.cpp src/testing.cpp src/PBC/test.cpp src/PBC/cuckoo.cpp src/PBC/mod.cpp src/Multi-PIR/pir.cpp src/Multi-PIR/test.cpp src/Multi-PIR/multipir.cpp src/OCC/occ.cpp src/OCC/test.cpp src/OCC/rbmatrices.cpp src/OCD/okvs.cpp src/OCD/test.cpp src/OCD/3H_GCT.cpp src/Parallel/opgpu.cu src/Parallel/parallel.cpp src/benchmark/test.cpp  src/Batch-PIR/main.cpp)
set(HEADER_FILES include/spiral.h include/core.h include/constants.h include/poly.h include/util.h include/client.h include/testing.h include/cuckoo.h include/mod.h include/pir.h include/multipir.h include/occ.h include/rbmatrices.h include/okvs.h include/3H_GCT.h include/parallel.h include/opgpu.cuh)
add_executable(spiral ${SRC_FILES} ${HEADER_FILES})


# Must specify:
# -DCMAKE_TOOLCHAIN_FILE=/home/ubuntu/vcpkg/scripts/buildsystems/vcpkg.cmake

set(CMAKE_PREFIX_PATH 
  "~/vcpkg/installed/x64-linux/share"
  "/app/vcpkg/installed/x64-linux/share"
)

# 链接hexl
find_package(HEXL CONFIG REQUIRED)

# 找不到时请去安装libssl-dev
find_package(OpenSSL REQUIRED)

find_package(CUDAToolkit REQUIRED)
if (CUDAToolkit_FOUND)
    message(STATUS "CUDAToolkit: found")
    message(STATUS "CUDAToolkit_VERSION=${CUDAToolkit_VERSION}")
    message(STATUS "CMAKE_CUDA_COMPILER=${CMAKE_CUDA_COMPILER}")
    message(STATUS  "CUDAToolkit_VERSION_MAJOR=${CUDAToolkit_VERSION_MAJOR}" )
    message(STATUS  "CUDAToolkit_VERSION_MINOR=${CUDAToolkit_VERSION_MINOR}" )
    message(STATUS  "CUDAToolkit_VERSION_PATCH=${CUDAToolkit_VERSION_PATCH}" )
    message(STATUS  "CUDAToolkit_BIN_DIR=${CUDAToolkit_BIN_DIR}" )
    message(STATUS  "CUDAToolkit_INCLUDE_DIRS=${CUDAToolkit_INCLUDE_DIRS}" )
    message(STATUS  "CUDAToolkit_LIBRARIES=${CUDAToolkit_LIBRARIES}" )
    message(STATUS  "CUDAToolkit_LIBRARY_DIR=${CUDAToolkit_LIBRARY_DIR}" )
    message(STATUS  "CUDAToolkit_LIBRARY_ROOT=${CUDAToolkit_LIBRARY_ROOT}" )
    message(STATUS  "CUDAToolkit_TARGET_DIR=${CUDAToolkit_TARGET_DIR}" )
    message(STATUS  "CUDAToolkit_NVCC_EXECUTABLE=${CUDAToolkit_NVCC_EXECUTABLE}" )
    # 如果导入目标存在，查询目标属性
    if (TARGET CUDA::cudart)
        get_target_property(_cudart_includes CUDA::cudart INTERFACE_INCLUDE_DIRECTORIES)
        get_target_property(_cudart_imported CUDA::cudart IMPORTED_LOCATION)
        if (NOT _cudart_imported)
            get_target_property(_cudart_imported CUDA::cudart IMPORTED_LOCATION_RELEASE)
        endif()
        get_target_property(_cudart_links CUDA::cudart INTERFACE_LINK_LIBRARIES)

        # 头文件、动态库、静态库
        message(STATUS "CUDA::cudart INTERFACE_INCLUDE_DIRECTORIES: ${_cudart_includes}")
        message(STATUS "CUDA::cudart IMPORTED_LOCATION: ${_cudart_imported}")
        message(STATUS "CUDA::cudart INTERFACE_LINK_LIBRARIES: ${_cudart_links}")
    else()
        message(STATUS "CUDA::cudart target not available")
    endif()
else()
    message(FATAL_ERROR "CUDAToolkit not found. Install CUDA Toolkit or set CUDAToolkit_ROOT / CUDA_PATH.")
endif()

target_include_directories(spiral PUBLIC include)
target_include_directories(spiral PUBLIC ~/vcpkg/installed/)
target_include_directories(spiral PUBLIC /app/vcpkg/installed/)
target_link_libraries(spiral
  PRIVATE
    CUDA::cudart
    HEXL::hexl
    OpenSSL::SSL
    OpenSSL::Crypto
)

add_compile_definitions(
    TEXP=$(TEXP)
    TEXPRIGHT=$(TEXPRIGHT)
    TCONV=$(TCONV)
    TGSW=$(TGSW)
    QPBITS=$(QPBITS)
    PVALUE=$(PVALUE)
    QNUMFIRST=$(QNUMFIRST)
    QNUMREST=$(QNUMREST)
    OUTN=$(OUTN)
)

