
######################################################################
### CMAKE VERSION
######################################################################

include(cmake/target_sources_local.cmake)
include(ExternalProject)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
message(STATUS ${CMAKE_MODULE_PATH})

#set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
#set(CMAKE_CXX_EXTENSIONS OFF)


######################################################################
### DISABLING IN-SOURCE BUILD
### (DO NOT MODIFY)
######################################################################
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

# Generate pif code
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RELEASE")
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
if(    NOT cmake_build_type_tolower STREQUAL "debug"
   AND NOT cmake_build_type_tolower STREQUAL "release"
   AND NOT cmake_build_type_tolower STREQUAL "relwithdebinfo")
  message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\". Allowed values are Debug, Release, RelWithDebInfo (case-insensitive).")
endif()

######################################################################
# Version information
######################################################################
include(version.cmake)
set(CYTNX_VERSION
  ${CYTNX_VERSION_MAJOR}.${CYTNX_VERSION_MINOR}.${CYTNX_VERSION_PATCH}
)

set(CYTNX_VARIANT_INFO "")

######################################################################
### PROJECT
######################################################################
project(CYTNX VERSION ${CYTNX_VERSION} LANGUAGES CXX C )
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(cytnx STATIC)
set_property(TARGET cytnx PROPERTY C_VISIBILITY_PRESET hidden)
set_property(TARGET cytnx PROPERTY VISIBILITY_INLINES_HIDDEN ON)
target_include_directories(cytnx
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_include_directories(cytnx SYSTEM
    PUBLIC
				${Boost_INCLUDE_DIRS}
)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wformat=0 -Wno-c++11-narrowing -w")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wformat=0 -w -fsized-deallocation")


target_compile_features(cytnx PUBLIC cxx_std_17)
#set_property(TARGET cytnx PROPERTY CXX_STANDARD 11)

add_subdirectory(src)

include(GNUInstallDirs)

######################################################################
### Find BLAS and LAPACK
######################################################################
set(BLA_VENDOR OpenBLAS)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(LAPACKE REQUIRED)
target_link_libraries(cytnx PUBLIC ${LAPACK_LIBRARIES})
target_link_libraries(cytnx PUBLIC ${LAPACKE_LIBRARIES})
message( STATUS "LAPACK found: ${LAPACK_LIBRARIES}" )
message( STATUS "LAPACKE found: ${LAPACKE_LIBRARIES}" )

find_package(Boost REQUIRED)
target_link_libraries(cytnx PUBLIC Boost::boost)

target_include_directories(cytnx PUBLIC /opt/homebrew/opt/openblas/include)
