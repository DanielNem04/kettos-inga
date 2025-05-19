cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

project(DoublePendulum LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(double_pendulum
    src/core.cpp
    src/simulation.cpp
    src/analysis.cpp
)

target_include_directories(double_pendulum PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

add_executable(pendulum_app app/main.cpp)

target_include_directories(pendulum_app PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(pendulum_app PRIVATE
    double_pendulum
)

install(TARGETS double_pendulum EXPORT DoublePendulumTargets
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(FILES
    include/double_pendulum/analysis.h
    include/double_pendulum/core.h
    include/double_pendulum/simulation.h
    DESTINATION include/double_pendulum
)

install(TARGETS pendulum_app
    RUNTIME DESTINATION bin
)

install(EXPORT DoublePendulumTargets
    FILE DoublePendulumConfig.cmake
    NAMESPACE DoublePendulum::
    DESTINATION share/DoublePendulum
)
