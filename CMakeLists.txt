project(qc)
cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_FLAGS "-std=c++0x")
aux_source_directory(./src SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})
find_package( OpenCV REQUIRED )
target_link_libraries( qc ${OpenCV_LIBS} )

