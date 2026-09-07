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
add_ig_test(test_perf_caches)

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

add_executable(diag_handle_order diag_handle_order.cpp)
target_link_libraries(diag_handle_order PRIVATE intersection_gen_lib)
target_include_directories(diag_handle_order PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(diag_all_violations diag_all_violations.cpp)
target_link_libraries(diag_all_violations PRIVATE intersection_gen_lib)
target_include_directories(diag_all_violations PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(diag_uturn_leads diag_uturn_leads.cpp)
target_link_libraries(diag_uturn_leads PRIVATE intersection_gen_lib)
target_include_directories(diag_uturn_leads PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(diag_family_ladder diag_family_ladder.cpp)
target_link_libraries(diag_family_ladder PRIVATE intersection_gen_lib)
target_include_directories(diag_family_ladder PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(diag_pair_cross diag_pair_cross.cpp)
target_link_libraries(diag_pair_cross PRIVATE intersection_gen_lib)
target_include_directories(diag_pair_cross PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(diag_fanout_order diag_fanout_order.cpp)
target_link_libraries(diag_fanout_order PRIVATE intersection_gen_lib)
target_include_directories(diag_fanout_order PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(diag_fan_natural diag_fan_natural.cpp)
target_link_libraries(diag_fan_natural PRIVATE intersection_gen_lib)
target_include_directories(diag_fan_natural PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(diag_fan_shape_space diag_fan_shape_space.cpp)
target_link_libraries(diag_fan_shape_space PRIVATE intersection_gen_lib)
target_include_directories(diag_fan_shape_space PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(diag_gen_time diag_gen_time.cpp)
target_link_libraries(diag_gen_time PRIVATE intersection_gen_lib)
target_include_directories(diag_gen_time PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(diag_uturn_nesting diag_uturn_nesting.cpp)
target_link_libraries(diag_uturn_nesting PRIVATE intersection_gen_lib)
target_include_directories(diag_uturn_nesting PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(diag_uturn_candidate_cross diag_uturn_candidate_cross.cpp)
target_link_libraries(diag_uturn_candidate_cross PRIVATE intersection_gen_lib)
target_include_directories(diag_uturn_candidate_cross PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(diag_single_cubic diag_single_cubic.cpp)
target_link_libraries(diag_single_cubic PRIVATE intersection_gen_lib)
target_include_directories(diag_single_cubic PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(diag_turning_excess diag_turning_excess.cpp)
target_link_libraries(diag_turning_excess PRIVATE intersection_gen_lib)
target_include_directories(diag_turning_excess PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(diag_fence diag_fence.cpp)
target_link_libraries(diag_fence PRIVATE intersection_gen_lib)
target_include_directories(diag_fence PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(diag_group_dir diag_group_dir.cpp)
target_link_libraries(diag_group_dir PRIVATE intersection_gen_lib)
target_include_directories(diag_group_dir PRIVATE ${CMAKE_SOURCE_DIR}/src)
