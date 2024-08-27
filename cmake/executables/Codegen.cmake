# Code generator.

add_executable(sc_make main/codegen/sc_make.cpp)
target_link_libraries(sc_make PRIVATE qontra)
target_include_directories(sc_make PRIVATE "include")
