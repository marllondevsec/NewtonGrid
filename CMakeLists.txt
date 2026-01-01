cmake_minimum_required(VERSION 3.15)
project(${PROJECT_NAME} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# targets
add_library(gravfield
    src/physics/mass.cpp
    src/physics/gravitational_field.cpp
    src/data/grid.cpp
    src/data/field_data.cpp
    src/io/export_csv.cpp
    src/io/export_vti.cpp
    src/utils/config.cpp
    src/utils/logger.cpp
    src/utils/math_utils.cpp
)

add_executable(gravfield_app src/main.cpp src/app/experiment_runner.cpp)
target_link_libraries(gravfield_app PRIVATE gravfield)

# tests (optional)
enable_testing()
add_subdirectory(src/tests)
