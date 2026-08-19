# Agent Rules

- 默认使用中文交互
- 只允许操作当前项目文件夹
- 每次修改完成后需按时间记录更新内容到docs/WORK_LOG.md
- Any code change that affects curve generation, constraints, ordering, shape preservation, or validation must be verified against all relevant constraints before it is considered complete.
- A fix is not complete if it only satisfies one failing scenario while breaking shape constraints, same-cluster non-endpoint non-intersection, fixed-shape preservation, Boundary avoidance, fence avoidance, obstacle avoidance, G1 continuity, or documented turn-shape requirements.
- For multi-constraint fixes, add or update regression tests that cover the original failure and the most likely interacting constraints. Run those tests before reporting completion.
- If the full test suite is too slow to run, run the focused regression tests that jointly cover the changed behavior and explicitly report any broader tests that were not run.
- Curve generation changes must preserve interactive performance: a single intersection should normally complete at millisecond scale, and the worst acceptable single-intersection generation time is under 15 seconds. Avoid exhaustive/global searches that can exceed this bound; when touching generation, constraints, ordering, shape preservation, or validation, verify representative single-intersection runtime before reporting completion.
- When adding/modifying/deleting functions and functional modules, comments and documentation need to be synchronized.
- If a conversation lasts longer than 20 minutes, it should automatically summarize and conclude, and provide suggestions for follow-up processing to prevent the process from going astray.
