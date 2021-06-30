cmake_minimum_required(VERSION 3.0)
project(mc-pricer CXX)

# Prefix
if (NOT PREFIX)
    set(PREFIX ${CMAKE_CURRENT_BINARY_DIR})
endif (NOT PREFIX)
set(CMAKE_INSTALL_PREFIX ${PREFIX} CACHE STRING "Installation prefix." FORCE)

if (CXX)
    set(CMAKE_CXX_COMPILER ${CXX})
endif (CXX)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# Release or Debug
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Debug' as none was specified.")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
endif ()


# Detect PNL
#find_package(Pnl REQUIRED)
#set(LIBS ${LIBS} ${PNL_LIBRARIES})
#include_directories(${PNL_INCLUDE_DIRS})
#if(${CMAKE_BUILD_TYPE} STREQUAL "Release")
 #   add_definitions(-DPNL_RANGE_CHECK_OFF)
#endif()

# Add jlparser


# Test the parser
add_executable(test_heston PriceCall_TG.cpp ModelHeston.cpp RandomNumberGenerator.cpp PathSimulator2F.cpp MCPricer.cpp Payoff.cpp  Utilities.cpp)
#target_link_libraries( test_heston ${LIBS})

#add_executable(mc-example src/mc_example.cpp src/PricingResults.cpp src/PnlVectToJson.cpp)
#target_link_libraries(mc-example jl_parser ${LIBS})

#add_executable(testAssetBPPast src/BlackScholesModel.cpp src/TestBSPast.cpp)
#target_link_libraries(testAssetBPPast ${LIBS})

#add_executable(testMontePast src/BlackScholesModel.cpp src/TestMonteCarloPasse.cpp src/Basket.cpp src/MonteCarlo.cpp)
#target_link_libraries(testMontePast ${LIBS})

# Print compiler flags
get_directory_property(DirDefs COMPILE_DEFINITIONS)
message(STATUS "COMPILE_DEFINITIONS = ${DirDefs}")
if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
    message(STATUS "C++ Flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_LINKER_EXE_FLAGS}" )
else()
    message(STATUS "C++ Flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_LINKER_EXE_FLAGS}" )
endif()
