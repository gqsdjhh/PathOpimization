cmake_minimum_required(VERSION 3.10)
project(a_star_grid)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)

find_package(Eigen3 REQUIRED)
find_package(osqp REQUIRED)
find_package(OsqpEigen REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE PYTHON_NUMPY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

add_executable(a_star_grid
    src/main.cpp
    src/grid.cpp
    src/astar.cpp
    src/visualize.cpp
    src/adaptive_segmentation.cpp
    src/polygon_corridor.cpp
    src/optimizer.cpp
)

target_include_directories(a_star_grid PRIVATE 
    ${PYTHON_NUMPY_INCLUDE_DIR}  # <--- 这里加上 NumPy 的路径
    ${Python3_INCLUDE_DIRS}      # 显式加上 Python 头文件路径 (保险起见)
)

target_link_libraries(a_star_grid 
Eigen3::Eigen 
osqp::osqp
OsqpEigen::OsqpEigen
Python3::Python
)
