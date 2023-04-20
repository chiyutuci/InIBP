# set library directories
set(POLY_SRC_DIR src)

# set library sources
aux_source_directory(${POLY_SRC_DIR} POLY_SRC)
add_library(poly ${POLY_SRC})
target_include_directories(poly PUBLIC ${POLY_SRC_DIR})

# add third party libraries
# flint library
find_library(FLINT_LIBRARY NAMES flint)
if (NOT FLINT_LIBRARY)
    message(FATAL_ERROR "flint library not found")
endif()
target_link_libraries(inibp ${FLINT_LIBRARY})