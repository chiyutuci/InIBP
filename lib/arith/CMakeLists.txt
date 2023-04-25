# set library sources
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} ARITH_SRC)
add_library(arith ${ARITH_SRC})
target_include_directories(arith PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# add third party libraries
# flint library
find_library(FLINT_LIBRARY NAMES flint)
if (NOT FLINT_LIBRARY)
    message(FATAL_ERROR "flint library not found")
endif()
target_link_libraries(inibp ${FLINT_LIBRARY})