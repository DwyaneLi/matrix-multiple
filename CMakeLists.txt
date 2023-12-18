cmake_minimum_required(VERSION 3.12)
project(Matrix_Multiple)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 添加可执行文件
add_executable(matrix_multiple matrix_multiple.cc)

# 查找并安装 OpenMP
find_package(OpenMP REQUIRED)

# 添加编译选项和链接库
 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
 set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

# 添加可执行文件依赖的头文件
target_include_directories(matrix_multiple PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(matrix_multiple PRIVATE -mavx512f)
target_link_libraries(matrix_multiple PRIVATE OpenMP::OpenMP_CXX)

install(TARGETS matrix_multiple DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin)
