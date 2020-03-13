cmake_minimum_required(VERSION 2.8.3)
project(cassie_estimation)

add_compile_options(-std=c++11)

## Find catkin macros and libraries
find_package(catkin REQUIRED COMPONENTS
  roscpp
  cassie_description
)

###################################
## catkin specific configuration ##
###################################
catkin_package(
  INCLUDE_DIRS include
  LIBRARIES ${PROJECT_NAME}
  CATKIN_DEPENDS roscpp cassie_description
  DEPENDS
)

###########
## Build ##
###########
add_definitions(-Wno-dev -Wno-deprecated-declarations)
include(common_flags)

## Specify additional locations of header files
## Your package locations should be listed before other locations
include_directories(
    include
    ${catkin_INCLUDE_DIRS}
)

option(USE_VECTORIZATION "Use vectorization" OFF)
if(USE_VECTORIZATION)
  common_add_vectorization_flags()
endif()

## Declare a C++ library
add_library(${PROJECT_NAME}
    src/contact_classifier.cpp
    src/heelspring_solver.cpp
    src/kinematics_hip_velocity_estimator.cpp
)
target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES})
