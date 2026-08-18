cmake_minimum_required(VERSION 3.15)
project(intersection_shape_generator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
option(INTERSECTION_BUILD_TESTS "Build project tests" ON)

# ── 性能优化: Debug模式使用-Og而非默认-O0 ──────────────────────
# -Og 在保留调试体验的同时进行轻量优化,使Debug性能从分钟级降至秒级。
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-g -Og" CACHE STRING "Debug flags" FORCE)

add_compile_definitions(PROJECT_ROOT_DIR="${PROJECT_SOURCE_DIR}")

set(DEPENDS_EXISTS "FALSE")
set(DEPENDS_DIR "${CMAKE_SOURCE_DIR}/depends")
if(EXISTS "${DEPENDS_DIR}")
    set(DEPENDS_EXISTS "TRUE")
    message(STATUS "Found [DEPENDS_DIR=${DEPENDS_DIR}], using it directly.")
    include_directories(${DEPENDS_DIR}/eigen/include/eigen3)
    link_directories(${DEPENDS_DIR}/eigen/lib)
else()
    set(DEPENDS_EXISTS "FALSE")
    message(STATUS "Not found [DEPENDS_DIR=${DEPENDS_DIR}], using <FetchContent> fetch libriaries.")
    include(FetchContent)

    ## eigen
    message(STATUS "Start fetch libriary: eigen ...")
    set(BUILD_TESTING OFF CACHE BOOL "Disable dependency tests" FORCE)
    option(EIGEN_BUILD_DOC "Enable creation of Eigen documentation" OFF)
    option(EIGEN_BUILD_PKGCONFIG "Build pkg-config .pc file for Eigen" OFF)
    option(EIGEN_SPLIT_LARGE_TESTS "Split large tests into smaller executables" OFF)
    FetchContent_Declare(eigen
            GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
            GIT_TAG        3.4.0
            GIT_SHALLOW    TRUE)
    FetchContent_MakeAvailable(eigen)
    if(INTERSECTION_BUILD_TESTS)
        set(BUILD_TESTING ON CACHE BOOL "Enable project tests" FORCE)
    endif()

endif()


add_subdirectory(src)
if(INTERSECTION_BUILD_TESTS)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests")
        enable_testing()
        add_subdirectory(tests)
    endif()
else()
    message(STATUS "Tests disabled by INTERSECTION_BUILD_TESTS=OFF")
endif()
