# Author: Suhas Vittal

if (NOT CPLEX_HOME)
    message(FATAL_ERROR "CPLEX is required to compile a library/application, 
    but CPLEX_HOME has not been set by the user.")
endif()

if (APPLE)
    set(ARCH_STRING "arm64_osx")
endif()

set(CPLEX_LIBRARY_DIR "${CPLEX_HOME}/cplex/lib/${ARCH_STRING}/static_pic")
set(CPLEX_INCLUDE_DIR "${CPLEX_HOME}/cplex/include/ilcplex")

set(CPLEX_LIBRARIES ${CPLEX_LIBRARY_DIR}/libcplex.a
                    ${CPLEX_LIBRARY_DIR}/libilocplex.a)
