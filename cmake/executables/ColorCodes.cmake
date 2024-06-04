# Color Codes

# PLACC independent executables:
add_executable(c2capacity main/color_codes/capacity.cpp)
add_executable(c2db main/color_codes/debugger.cpp)

target_link_libraries(c2capacity PRIVATE qontra)
target_link_libraries(c2db PRIVATE qontra)

# PLACC dependent executables:
set(PLACC_FILES
        src/placc/cx.cpp
        src/placc/fpn.cpp
        src/placc/fpn/metric.cpp
        src/placc/fpn/phase_one.cpp
        src/placc/fpn/phase_two.cpp
        src/placc/tree.cpp)
add_library(placc ${PLACC_FILES})
target_link_libraries(placc PUBLIC qontra)

add_executable(c2fpnmake main/color_codes/fpnmake.cpp)
target_link_libraries(c2fpnmake PRIVATE placc)
