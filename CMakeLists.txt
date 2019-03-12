
project(reconstruct)
cmake_minimum_required(VERSION 3.10)
find_package(CGAL)
find_package( Boost 1.65.1 COMPONENTS program_options REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )
add_executable(reconstruct main.cpp)
target_link_libraries(reconstruct CGAL::CGAL ${Boost_LIBRARIES})

find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater) 
if (EIGEN3_FOUND) 
   include( ${EIGEN3_USE_FILE} ) 
endif() 