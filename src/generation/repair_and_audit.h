#pragma once

// 修复与审计模块：集中修复预算、影响闭包、事务提交和最终状态投影。
#include "generation/atomic_curve_patch.h"
#include "generation/final_curve_auditor.h"
#include "generation/final_pair_auditor.h"
#include "generation/physical_repair_coordinator.h"
#include "generation/repair_budget.h"
#include "generation/repair_impact_closure.h"

