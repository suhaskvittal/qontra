# Color Codes

# PLACC independent executables:
add_executable(c2capacity main/color_codes/capacity.cpp)
add_executable(c2db main/color_codes/debugger.cpp)
add_executable(mat2 main/color_codes/mat2.cpp)

target_link_libraries(c2capacity PRIVATE qontra)
target_link_libraries(c2db PRIVATE qontra)
target_link_libraries(mat2 PRIVATE qontra)

add_executable(c2make main/color_codes/make_code.cpp
                        src/codegen/tiling.cpp
                        src/codegen/convert.cpp
                        src/codegen/monte_carlo.cpp)
target_link_libraries(c2make PRIVATE qontra)
target_include_directories(c2make PRIVATE "include")

