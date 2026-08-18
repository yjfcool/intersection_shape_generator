set(ISG_PUBLIC_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

function(isg_configure_target target)
    target_include_directories(${target} PUBLIC ${ISG_PUBLIC_INCLUDE_DIR})
    if (NOT DEPENDS_EXISTS STREQUAL "TRUE")
        target_link_libraries(${target} PUBLIC Eigen3::Eigen)
    endif ()
    target_compile_options(${target} PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra>
        $<$<CXX_COMPILER_ID:MSVC>:/W3>
        $<$<CONFIG:Debug>:-g -Og -DNDEBUG>
        $<$<CONFIG:Release>:-O3 -DNDEBUG>
        $<$<CONFIG:RelWithDebInfo>:-O2 -g -DNDEBUG>)
endfunction()

add_library(isg_core STATIC

        domain/scene_view.cpp
        geometry/predicates.cpp
        utils/clipper.cpp
        curve/bezier.cpp
        curve/hermite_init.cpp
        curve/curve_utils.cpp
        optimizer/sdf_field.cpp
        optimizer/penalty_cost.cpp
        optimizer/lbfgs_solver.cpp
        constraints/infeasibility_detector.cpp
        constraints/cluster_order.cpp
        constraints/constraint_evaluator.cpp
        constraints/shape_constraint.cpp
        constraints/cluster_constraint.cpp
        constraints/uturn_envelope_constraint.cpp
        toolkits/toolkits.cpp
        preprocessing/uturn_family_builder.cpp
    preprocessing/crosswalk_clearance_calculator.cpp
    initialization/curve_initializers.cpp
    ordering/generation_planner.cpp
    generation/curve_sampling.cpp
    generation/uturn_multi_constraint_search.cpp
    generation/uturn_shape.cpp
    generation/segmented_uturn_candidate_search.cpp
    generation/connectivity_generation_context.cpp
    generation/candidate_generation.cpp
    generation/curve_processing.cpp
    generation/repair_and_audit.cpp
    generator/connectivity_generator.cpp
    generator/connectivity_generation_session.cpp
        generator/polygon_builder.cpp
        generator/edge_line_generator.cpp
        io/iodata_json.cpp
        io/iodata_shapefile.cpp)
isg_configure_target(isg_core)


# 保留原公共目标名，作为稳定的 API 门面。
add_library(intersection_gen_lib STATIC intersection_shape_generator.cpp)
isg_configure_target(intersection_gen_lib)
target_link_libraries(intersection_gen_lib PUBLIC isg_core)

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")
    add_executable(intersection_gen main.cpp)
    target_link_libraries(intersection_gen PRIVATE intersection_gen_lib)
    target_compile_options(intersection_gen PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra>
        $<$<CXX_COMPILER_ID:MSVC>:/W3>)
endif ()
