cmake_minimum_required(VERSION 3.16)

# set c++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# add compile options

# project name
project(inibp VERSION 0.1.0)

# set project directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(INIBP_SRC_DIR src)
set(INIBP_LIB_DIR lib)
set(INIBP_INCLUDE_DIR include)
set(INIBP_LIBRARIES arith)

# add source files
aux_source_directory(${INIBP_SRC_DIR} INIBP_SRC)
add_executable(inibp ${INIBP_SRC})
target_include_directories(inibp PUBLIC ${INIBP_SRC_DIR} ${INIBP_INCLUDE_DIR} ${INIBP_LIB_DIR})

# add libraries
foreach (LIBRARY ${INIBP_LIBRARIES})
    add_subdirectory(${INIBP_LIB_DIR}/${LIBRARY})
endforeach ()
target_link_libraries(inibp ${INIBP_LIBRARIES})

# add third party libraries
# yaml-cpp library
find_library(YAML_LIBRARY NAMES yaml-cpp)
if (NOT YAML_LIBRARY)
    message(FATAL_ERROR "yaml-cpp library not found")
endif ()
target_link_libraries(inibp ${YAML_LIBRARY})

#ginac library
find_library(GINAC_LIBRARY NAMES ginac)
if (NOT GINAC_LIBRARY)
    message(FATAL_ERROR "ginac library not found")
endif ()
target_link_libraries(inibp ${GINAC_LIBRARY})