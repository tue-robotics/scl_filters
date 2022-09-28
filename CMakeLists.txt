cmake_minimum_required(VERSION 2.8.3)
project(scl_filters)

add_compile_options(-Wall -Werror=all)
add_compile_options(-Wextra -Werror=extra)

find_package(catkin REQUIRED)

catkin_package(
 INCLUDE_DIRS include
 LIBRARIES scl_filters
)

include_directories(include)

add_library(scl_filters 	 	 src/DFirstOrderLowpass.cpp
								 src/DPD.cpp
                                 src/DSecondOrderLowpass.cpp
                                 src/DWeakIntegrator.cpp
                                 src/DLeadLag.cpp
                                 src/DPID.cpp
                                 src/DSkewedNotch.cpp
                                 src/Polynomial.cpp
                                 )
