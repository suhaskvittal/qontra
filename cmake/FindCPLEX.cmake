# Author: Suhas Vittal

if (NOT CPLEX_HOME)
    set(CPLEX_HOME "/Applications/CPLEX_Studio2211") # for MacOS
endif()

if (APPLE)
    set(ARCH_STRING "arm64_osx")
endif()

set(CPLEX_LIB_PREFIX "${CPLEX_HOME}/cplex/lib/${ARCH_STRING}/static_pic")
set(CPLEX_INCLUDE_DIR "${CPLEX_HOME}/cplex/include/ilcplex")
