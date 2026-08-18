#include "intersection_shape_generator.h"
#include "constraints/cluster_constraint.h"
#include "constraints/constraint_evaluator.h"
#include "constraints/shape_constraint.h"
#include "constraints/uturn_envelope_constraint.h"
#include "domain/scene_context.h"
#include "generation/connectivity_generation_context.h"
#include "generation/candidate_selector.h"
#include "generation/bounded_uturn_candidate_search.h"
#include "generation/curve_post_processor.h"
#include "generation/curve_optimizer.h"
#include "generation/optimization_result_processor.h"
#include "generation/physical_repair_coordinator.h"
#include "generation/final_curve_auditor.h"
#include "generation/repair_budget.h"
#include "generation/repair_impact_closure.h"
#include "generation/atomic_curve_patch.h"
#include "generation/final_pair_auditor.h"
#include "initialization/curve_initializer_registry.h"
#include "preprocessing/uturn_family_builder.h"
#include "preprocessing/crosswalk_clearance_calculator.h"
#include "toolkits/toolkits.h"
#include "geometry/predicates.h"
#include "ordering/generation_planner.h"
#include "generator/edge_line_generator.h"
#include "generator/polygon_builder.h"
#include "io/iodata_json.h"
#include "io/iodata_shapefile.h"

int main() {
    isg::IntersectionShapeGenerator generator;
    isg::EdgeLineGenerator edge_generator;
    isg::IntersectionAreaBuilder area_builder(0.5, 4.0);
    (void)generator;
    (void)edge_generator;
    (void)area_builder;
    return 0;
}
