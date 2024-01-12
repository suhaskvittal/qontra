# Author:   Suhas Vittal
# date:     25 December 2023

add_executable(converter main/converter.cpp)
add_executable(syndromes main/syndromes.cpp)

target_link_libraries(converter PRIVATE qontra)
target_link_libraries(syndromes PRIVATE qontra)
