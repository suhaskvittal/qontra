# Color Codes

add_executable(c2capacity main/color_codes/capacity.cpp)
add_executable(c2db main/color_codes/debugger.cpp)

target_link_libraries(c2capacity PRIVATE qontra)
target_link_libraries(c2db PRIVATE qontra)

