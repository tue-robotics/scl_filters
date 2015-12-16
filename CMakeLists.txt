cmake_minimum_required(VERSION 2.8.3)
project(scl_filters)
find_package(catkin REQUIRED)

catkin_package(
    INCLUDE_DIRS include
    LIBRARIES scl_filters
)

include_directories(include)

add_library(scl_filters
    include/scl/polynomial/Polynomial.hpp
    include/scl/filters/DSecondOrderLowpass.hpp
    include/scl/filters/DPID.hpp
    include/scl/filters/DPD.hpp
    include/scl/filters/DSkewedNotch.hpp
    include/scl/filters/DFirstOrderLowpass.hpp
    include/scl/filters/DiscretizationMethods.hpp
    include/scl/filters/DWeakIntegrator.hpp
    include/scl/filters/DLeadLag.hpp

    src/DFirstOrderLowpass.cpp
    src/DPD.cpp
    src/DSecondOrderLowpass.cpp
    src/DWeakIntegrator.cpp
    src/DLeadLag.cpp
    src/DPID.cpp
    src/DSkewedNotch.cpp
    src/Polynomial.cpp
)
