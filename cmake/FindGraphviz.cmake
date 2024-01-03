# Author: Suhas Vittal

if (NOT GRAPHVIZ_INSTALL_DIR)
    message(FATAL_ERROR "Graphviz is required to compile a library/application, 
    but GRAPHVIZ_INSTALL_DIR has not been set by the user.")
endif()

set(GRAPHVIZ_INCLUDE_DIR "${GRAPHVIZ_INSTALL_DIR}/include")
set(GRAPHVIZ_LIBRARY_DIR "${GRAPHVIZ_INSTALL_DIR}/lib")

set(GRAPHVIZ_LIBRARIES ${GRAPHVIZ_LIBRARY_DIR}/libcgraph.dylib
                        ${GRAPHVIZ_LIBRARY_DIR}/libgvc.dylib)
