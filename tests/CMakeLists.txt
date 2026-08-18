include(FetchContent)

# Catch2 v3 (header+source, MIT)
FetchContent_Declare(Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.5.4
    GIT_SHALLOW    TRUE)
FetchContent_MakeAvailable(Catch2)

# Helper to register a test executable
function(add_ig_test name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PRIVATE
        intersection_gen_lib
        Catch2::Catch2WithMain
        )
    target_include_directories(${name} PRIVATE
        ${CMAKE_SOURCE_DIR}/src)
    add_test(NAME ${name} COMMAND ${name})
endfunction()

enable_testing()

add_ig_test(test_bezier)
add_ig_test(test_sdf)
add_ig_test(test_optimizer)
add_ig_test(test_cluster)
add_ig_test(test_edge_line)
add_ig_test(test_intersection_input_regression)
add_ig_test(test_infeasibility)
add_ig_test(test_integration)
add_ig_test(test_scenarios)
add_ig_test(test_refactored_modules)

# Compile-only public-header self-containment check; no Catch2 dependency needed.
add_executable(test_public_headers test_public_headers.cpp)
target_link_libraries(test_public_headers PRIVATE intersection_gen_lib)
target_include_directories(test_public_headers PRIVATE ${CMAKE_SOURCE_DIR}/src)
add_test(NAME test_public_headers COMMAND test_public_headers)

# Standalone diagnostic executable (not a test, just runs and prints)
add_executable(diag_uturn diag_uturn.cpp)
target_link_libraries(diag_uturn PRIVATE intersection_gen_lib)
target_include_directories(diag_uturn PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(diag_uturn_geom diag_uturn_geom.cpp)
target_link_libraries(diag_uturn_geom PRIVATE intersection_gen_lib)
target_include_directories(diag_uturn_geom PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(diag_cluster_pairs diag_cluster_pairs.cpp)
target_link_libraries(diag_cluster_pairs PRIVATE intersection_gen_lib)
target_include_directories(diag_cluster_pairs PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(diag_pair_isect diag_pair_isect.cpp)
target_link_libraries(diag_pair_isect PRIVATE intersection_gen_lib)
target_include_directories(diag_pair_isect PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(diag_focus diag_focus.cpp)
target_link_libraries(diag_focus PRIVATE intersection_gen_lib)
target_include_directories(diag_focus PRIVATE ${CMAKE_SOURCE_DIR}/src)
