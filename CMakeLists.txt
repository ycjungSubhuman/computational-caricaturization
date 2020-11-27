project(computational-caricaturization C CXX)
cmake_minimum_required(VERSION 3.11)

find_package(Eigen3 REQUIRED)
find_package(Threads REQUIRED)

add_library(compcari src/compcari.cpp)
target_link_libraries(compcari PUBLIC Eigen3::Eigen Threads::Threads)
target_include_directories(compcari PRIVATE 3rdparty/libigl/include include)

add_executable(main demo/main.cpp)
target_include_directories(main PRIVATE 3rdparty/libigl/include include)
target_link_libraries(main PRIVATE compcari)