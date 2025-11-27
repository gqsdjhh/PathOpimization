cmake_minimum_required(VERSION 3.10)
project(a_star_grid)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)

find_package(Eigen3 REQUIRED)
find_package(osqp REQUIRED)
find_package(OsqpEigen REQUIRED)

add_executable(a_star_grid
    src/main.cpp
    src/grid.cpp
    src/astar.cpp
    src/visualize.cpp
    src/adaptive_segmentation.cpp
    src/polygon_corridor.cpp
    src/optimizer.cpp
)

target_link_libraries(a_star_grid 
Eigen3::Eigen 
osqp::osqp
OsqpEigen::OsqpEigen
)
