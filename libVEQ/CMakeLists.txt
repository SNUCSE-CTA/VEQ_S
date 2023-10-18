add_library(libVEQ STATIC)

target_sources(libVEQ PRIVATE
    libVEQ.cpp
    libVEQ.hpp
)

target_include_directories(libVEQ PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)
