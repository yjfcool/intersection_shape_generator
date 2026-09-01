#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""统一测试报表：把四类测试结果按「检查类别」汇总成 Excel（每类别一个 sheet）。

覆盖的四类测试（与 docs/WORK_LOG.md 的验证小节一一对应）：
  全量审计  tests/diag_all_violations —— 逐份数据集生成后复算全部约束，输出 findings CSV
  回归套件  tests/test_*             —— Catch2 套件，用 XML reporter 输出后解析
  焦点用例  tools/focus_cases.txt    —— 高价值用例逐条单跑，各自一份 XML，便于定位断言
  性能      tests/diag_gen_time      —— 单路口生成耗时（15s 硬上限），输出 findings CSV

所有结果先归一到同一张发现表（列见 FIELDS），再由 tools/xlsx_writer.py 写出：
  汇总 / 按类别统计 / 按文件统计 / 类别×文件矩阵 / 各类别明细 / 全部明细 / 用例结果
  以及可选的「基线对比」（--baseline 指向上一轮的 findings.csv，标注 新增/保持/已消除）。

典型用法：
  python3 tools/test_report.py                          # 跑全部四类并出报表
  python3 tools/test_report.py --stages audit,perf      # 只跑审计与性能
  python3 tools/test_report.py --skip-run               # 全部复用已有产物，只重出报表
  python3 tools/test_report.py --skip-run audit,perf    # 只复用审计与性能，回归/焦点实跑
  python3 tools/test_report.py --baseline reports/findings_base.csv
"""

import argparse
import csv
import os
import re
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import xlsx_writer as X  # noqa: E402  （同目录模块，需先补 sys.path）

# 发现表列定义：C++ 侧（diag_all_violations / diag_gen_time）写出的 CSV 必须与此一致。
FIELDS = ["source", "category", "check", "dataset", "subject", "severity",
          "metric", "value", "threshold", "detail", "location"]

SEVERITY_ORDER = {"error": 0, "warn": 1, "info": 2}

# 检查类别 id -> 中文名。sheet 顺序即此处顺序（按修复优先级从形态到辅助信息排列）。
CATEGORY_LABELS = OrderedDict([
    ("cluster_cross", "同簇相交"),
    ("cluster_precursor", "同簇相交前兆"),
    ("uturn_shape", "调头形态"),
    ("ordinary_shape", "常规形态"),
    ("crosswalk", "人行横道"),
    ("fixed_shape", "固有形态保持"),
    ("g1", "G1连续"),
    ("curvature", "曲率"),
    ("obstacle", "障碍避让"),
    ("fence", "围栏越界"),
    ("boundary", "边界避让"),
    ("road_edge", "路缘避让"),
    ("self_intersect", "自交"),
    ("area", "路口面"),
    ("lane_edge", "车道边线"),
    ("infeasibility", "不可行判定"),
    ("topology", "拓扑校验"),
    ("generation", "生成失败"),
    ("generator_report", "生成器自报"),
    ("perf", "性能"),
    ("suite", "套件与用例"),
    ("other", "其他"),
])

SOURCE_LABELS = OrderedDict([
    ("audit", "全量审计"),
    ("regression", "回归套件"),
    ("focus", "焦点用例"),
    ("perf", "性能"),
])

# Catch2 标签 -> 检查类别。按此顺序匹配，先命中者胜（越具体的标签排越前）。
TAG_CATEGORY = OrderedDict([
    ("merge-funnel", "cluster_cross"), ("cluster", "cluster_cross"),
    ("uturn", "uturn_shape"), ("crosswalk", "crosswalk"),
    ("fixed-geometry", "fixed_shape"), ("fixed_shape", "fixed_shape"),
    ("g1", "g1"), ("curvature", "curvature"),
    ("obstacle", "obstacle"), ("bypass", "obstacle"),
    ("fence", "fence"), ("boundary", "boundary"),
    ("road_edge", "road_edge"), ("roadedge", "road_edge"), ("clear", "road_edge"),
    ("area", "area"), ("edge_line", "lane_edge"), ("edge", "lane_edge"),
    ("shape", "ordinary_shape"), ("infeasibility", "infeasibility"),
    ("topo_block", "topology"), ("perf", "perf"),
])

# Catch2 套件（非 Catch2 的编译期自检单独跑，见 PLAIN_SUITES）
CATCH_SUITES = ["test_bezier", "test_sdf", "test_optimizer", "test_cluster",
                "test_edge_line", "test_intersection_input_regression",
                "test_infeasibility", "test_integration", "test_scenarios",
                "test_refactored_modules"]
PLAIN_SUITES = ["test_public_headers"]

MAX_FINDINGS_PER_CASE = 20  # 单条用例最多展开多少条失败断言，避免明细表被一条用例淹没

# 用例名关键词 -> 检查类别（标签给不出类别时的兜底，顺序即优先级）。
NAME_CATEGORY = [
    ("crosswalk", "crosswalk"), ("u-turn", "uturn_shape"), ("uturn", "uturn_shape"),
    ("curvature", "curvature"), ("cluster", "cluster_cross"),
    ("fixed", "fixed_shape"), ("fence", "fence"), ("obstacle", "obstacle"),
    ("boundary", "boundary"), ("road edge", "road_edge"), ("roadedge", "road_edge"),
    ("clearance", "road_edge"), ("self-intersect", "self_intersect"),
    ("g1", "g1"), ("continuity", "g1"), ("elapsed", "perf"),
    ("area", "area"), ("lane edge", "lane_edge"), ("edge line", "lane_edge"),
    ("arch", "ordinary_shape"), ("shape", "ordinary_shape"),
    ("infeasib", "infeasibility"), ("topolog", "topology"),
]

def label_of(category):
    """未登记的类别 id 原样显示，避免新增类别时静默归入“其他”。"""
    return CATEGORY_LABELS.get(category, category)


def new_finding(source, category, check, dataset, subject, severity,
                metric="", value="", threshold="", detail="", location=""):
    return {"source": source, "category": category, "check": check,
            "dataset": dataset or "-", "subject": subject or "-",
            "severity": severity, "metric": metric, "value": value,
            "threshold": threshold, "detail": detail, "location": location}


def read_findings_csv(path):
    """读取 C++ 工具写出的发现 CSV；缺列按空串补齐。"""
    rows = []
    if not os.path.isfile(path):
        return rows
    with open(path, "r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            rows.append({k: (row.get(k) or "") for k in FIELDS})
    return rows


def write_findings_csv(path, findings):
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS, extrasaction="ignore")
        writer.writeheader()
        for item in findings:
            writer.writerow(item)


def dataset_stems(datas_dir):
    """datas 目录下的数据集名（去掉 .json），用于把用例名/标签映射到具体文件。"""
    stems = []
    if os.path.isdir(datas_dir):
        for name in sorted(os.listdir(datas_dir)):
            if name.endswith(".json"):
                stems.append(name[:-5])
    return stems


def guess_dataset(text, tags, stems):
    """从用例名与标签里猜数据集：优先精确等于某个数据集名，其次前缀匹配最长者。"""
    tokens = [t for t in tags]
    tokens += re.findall(r"[0-9A-Za-z_][0-9A-Za-z_\-]*", text or "")
    for token in tokens:
        if token in stems:
            return token + ".json"
    best = ""
    for token in tokens:
        if not re.match(r"^\d{6,}", token):
            continue
        for stem in stems:
            if stem.startswith(token) and len(stem) > len(best):
                best = stem
    return best + ".json" if best else "-"


def guess_category(tags, name=""):
    """先按标签定类别，标签不足时再按用例名关键词兜底（否则大量用例只能归入“套件”）。"""
    for tag in tags:
        low = tag.lower()
        if low in TAG_CATEGORY:
            return TAG_CATEGORY[low]
    for tag in tags:
        low = tag.lower()
        for key, category in TAG_CATEGORY.items():
            if key in low:
                return category
    low_name = (name or "").lower()
    for key, category in NAME_CATEGORY:
        if key in low_name:
            return category
    return "suite"


def rel_path(path):
    """把绝对路径压成仓库内相对路径，便于“按文件”定位。"""
    if not path:
        return ""
    try:
        rel = os.path.relpath(path, ROOT)
    except ValueError:
        return path
    return path if rel.startswith("..") else rel


def run_cmd(cmd, log_path):
    """执行子进程，stdout/stderr 一并落盘，返回 (返回码, 耗时秒)。"""
    t0 = time.time()
    with open(log_path, "w", encoding="utf-8", errors="replace") as log:
        log.write("$ " + " ".join(cmd) + "\n\n")
        log.flush()
        try:
            code = subprocess.call(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=ROOT)
        except OSError as exc:
            log.write("\n[执行失败] %s\n" % exc)
            code = 127
    return code, time.time() - t0


def _loc(node):
    name = rel_path(node.get("filename"))
    line = node.get("line")
    return "%s:%s" % (name, line) if name and line else name


def _detail(sections, head, expanded):
    parts = []
    if sections:
        parts.append("[" + " > ".join(s for s in sections if s) + "]")
    if head:
        parts.append(head)
    if expanded and expanded != head:
        parts.append("=> " + expanded)
    return " ".join(parts)


def collect_results(node, sections=(), fails=None, skips=None):
    """递归收集用例内的失败断言与 SKIP，保留 SECTION 路径作为定位上下文。"""
    if fails is None:
        fails, skips = [], []
    for child in node:
        if child.tag == "Section":
            collect_results(child, sections + (child.get("name") or "",), fails, skips)
        elif child.tag == "Expression":
            if child.get("success") == "false":
                kind = re.sub(r"[^a-z0-9]+", "_",
                              (child.get("type") or "assert").lower())
                fails.append((kind, _loc(child),
                              _detail(sections, (child.findtext("Original") or "").strip(),
                                      (child.findtext("Expanded") or "").strip())))
            collect_results(child, sections, fails, skips)
        elif child.tag in ("Failure", "Exception", "FatalErrorCondition"):
            fails.append((child.tag.lower(), _loc(child),
                          _detail(sections, (child.text or "").strip(), "")))
        elif child.tag == "Skip":
            skips.append((_loc(child), _detail(sections, (child.text or "").strip(), "")))
    return fails, skips


def parse_catch2_xml(path, source, stems):
    """解析一份 Catch2 XML，返回 (findings, cases)。

    findings 只收失败/跳过（通过的用例不进明细，避免正常结果淹没问题）；
    cases 收全部用例结果，用于「用例结果」sheet 的通过率统计。
    """
    findings, cases = [], []
    if not os.path.isfile(path):
        return findings, cases
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as exc:
        findings.append(new_finding(source, "suite", "suite.xml_parse_error", "-",
                                    os.path.basename(path), "error",
                                    detail="XML 解析失败（进程可能崩溃）: %s" % exc,
                                    location=rel_path(path)))
        return findings, cases

    suite = root.get("name") or os.path.basename(path)
    for tc in root.iter("TestCase"):
        name = tc.get("name") or "-"
        tags = re.findall(r"\[([^\]]+)\]", tc.get("tags") or "")
        category = guess_category(tags, name)
        dataset = guess_dataset(name, tags, stems)
        location = _loc(tc)
        overall = tc.find("OverallResult")
        success = overall is not None and overall.get("success") == "true"
        skips = int(overall.get("skips") or 0) if overall is not None else 0
        duration = float(overall.get("durationInSeconds") or 0.0) \
            if overall is not None else 0.0
        fails, skip_msgs = collect_results(tc)
        result = "fail" if not success else ("skip" if skips else "pass")
        cases.append({"suite": suite, "source": source, "case": name,
                      "tags": " ".join("[%s]" % t for t in tags), "result": result,
                      "failed": len(fails), "skips": skips, "duration": duration,
                      "dataset": dataset, "category": category, "location": location})

        for i, (kind, loc, detail) in enumerate(fails):
            if i >= MAX_FINDINGS_PER_CASE:
                findings.append(new_finding(
                    source, category, "case.truncated", dataset, name, "info",
                    detail="该用例另有 %d 条失败断言未展开，详见 %s"
                           % (len(fails) - i, rel_path(path)), location=location))
                break
            findings.append(new_finding(source, category, "case." + kind, dataset,
                                        name, "error", detail=detail,
                                        location=loc or location))
        if not success and not fails:
            findings.append(new_finding(source, category, "case.failed", dataset, name,
                                        "error", detail="用例失败但 XML 无断言明细",
                                        location=location))
        for loc, detail in skip_msgs:
            findings.append(new_finding(source, category, "case.skipped", dataset, name,
                                        "warn", detail=detail or "用例被 SKIP",
                                        location=loc or location))
    return findings, cases


def binary_path(build_dir, name):
    return os.path.join(build_dir, "tests", name)


def reuse(args, stage):
    """该阶段是否复用已有产物（--skip-run 不带值=全部阶段，带值=指定子集）。"""
    spec = (args.skip_run or "").strip()
    if not spec:
        return False
    if spec == "all":
        return True
    return stage in [s.strip() for s in spec.split(",")]


def missing_binary(source, name):
    return new_finding(source, "suite", "suite.missing_binary", "-", name, "error",
                       detail="未找到可执行文件，请先构建：cmake --build <build> --target " + name,
                       location="tests/CMakeLists.txt")


def stage_audit(args, ctx):
    """全量审计：一次生成 + 全约束复算，CSV 由 C++ 侧直接产出。"""
    findings = []
    csv_path = os.path.join(ctx["out"], "findings_audit.csv")
    exe = binary_path(args.build_dir, "diag_all_violations")
    if not reuse(args, "audit"):
        if not os.path.isfile(exe):
            return [missing_binary("audit", "diag_all_violations")]
        cmd = [exe, args.datas, "--csv", csv_path]
        if args.only:
            cmd += ["--only", args.only]
        if args.road_edge_clearance is not None:
            cmd += ["--road-edge-clearance", str(args.road_edge_clearance)]
        code, secs = run_cmd(cmd, os.path.join(ctx["logs"], "audit.log"))
        print("  全量审计 完成 exit=%d %.1fs" % (code, secs), flush=True)
        if code > 1:  # 返回码 1 = 存在违约（正常业务结果），>1 才是工具异常
            findings.append(new_finding(
                "audit", "suite", "suite.tool_error", "-", "diag_all_violations",
                "error", detail="审计工具异常退出 exit=%d，详见 logs/audit.log" % code,
                location="tests/diag_all_violations.cpp"))
    findings += read_findings_csv(csv_path)
    return findings


def stage_perf(args, ctx):
    """性能：单路口生成耗时基准，与审计里的顺带计时相互印证。"""
    findings = []
    csv_path = os.path.join(ctx["out"], "findings_perf.csv")
    exe = binary_path(args.build_dir, "diag_gen_time")
    if not reuse(args, "perf"):
        if not os.path.isfile(exe):
            return [missing_binary("perf", "diag_gen_time")]
        cmd = [exe, "--dir", args.datas, "--csv", csv_path,
               "--repeat", str(args.repeat)]
        code, secs = run_cmd(cmd, os.path.join(ctx["logs"], "perf.log"))
        print("  性能 完成 exit=%d %.1fs" % (code, secs), flush=True)
        if code != 0:
            findings.append(new_finding(
                "perf", "suite", "suite.tool_error", "-", "diag_gen_time", "error",
                detail="性能工具异常退出 exit=%d，详见 logs/perf.log" % code,
                location="tests/diag_gen_time.cpp"))
    findings += read_findings_csv(csv_path)
    return findings


def stage_regression(args, ctx):
    """回归套件：逐个 Catch2 套件跑 XML reporter，另跑编译期自检可执行文件。"""
    findings, cases = [], []
    suites = args.suites.split(",") if args.suites else CATCH_SUITES
    for suite in [s.strip() for s in suites if s.strip()]:
        exe = binary_path(args.build_dir, suite)
        xml_path = os.path.join(ctx["xml"], suite + ".xml")
        if not reuse(args, "regression"):
            if not os.path.isfile(exe):
                findings.append(missing_binary("regression", suite))
                continue
            code, secs = run_cmd([exe, "--reporter", "xml", "--out", xml_path,
                                  "--durations", "yes"],
                                 os.path.join(ctx["logs"], suite + ".log"))
            print("  回归套件 %s exit=%d %.1fs" % (suite, code, secs), flush=True)
        f, c = parse_catch2_xml(xml_path, "regression", ctx["stems"])
        findings += f
        cases += c

    for suite in PLAIN_SUITES:
        exe = binary_path(args.build_dir, suite)
        if reuse(args, "regression"):
            continue
        if not os.path.isfile(exe):
            findings.append(missing_binary("regression", suite))
            continue
        code, secs = run_cmd([exe], os.path.join(ctx["logs"], suite + ".log"))
        print("  回归套件 %s exit=%d %.1fs" % (suite, code, secs), flush=True)
        cases.append({"suite": suite, "source": "regression", "case": suite,
                      "tags": "[compile-check]",
                      "result": "pass" if code == 0 else "fail", "failed": 0,
                      "skips": 0, "duration": secs, "dataset": "-",
                      "category": "suite", "location": "tests/%s.cpp" % suite})
        if code != 0:
            findings.append(new_finding("regression", "suite", "suite.plain_failed", "-",
                                        suite, "error",
                                        detail="非 Catch2 自检程序退出码 %d" % code,
                                        location="tests/%s.cpp" % suite))
    return findings, cases


def load_focus_list(path):
    """读取焦点用例清单，返回 [(可执行文件, 用例名)]。"""
    items = []
    if not os.path.isfile(path):
        return items
    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            suite, _, case = line.partition("|")
            if case.strip():
                items.append((suite.strip(), case.strip()))
    return items


def stage_focus(args, ctx):
    """焦点用例：逐条单跑，各自一份 XML，失败断言可精确到 SECTION。"""
    findings, cases = [], []
    items = load_focus_list(args.focus_list)
    if not items:
        findings.append(new_finding("focus", "suite", "suite.focus_list_missing", "-",
                                    rel_path(args.focus_list), "warn",
                                    detail="焦点用例清单为空或不存在，已跳过该阶段",
                                    location=rel_path(args.focus_list)))
        return findings, cases
    for index, (suite, case) in enumerate(items, 1):
        exe = binary_path(args.build_dir, suite)
        xml_path = os.path.join(ctx["xml"], "focus_%02d.xml" % index)
        if not reuse(args, "focus"):
            if not os.path.isfile(exe):
                findings.append(missing_binary("focus", suite))
                continue
            code, secs = run_cmd([exe, case, "--reporter", "xml", "--out", xml_path,
                                  "--durations", "yes"],
                                 os.path.join(ctx["logs"], "focus_%02d.log" % index))
            print("  焦点用例 %d/%d exit=%d %.1fs %s"
                  % (index, len(items), code, secs, case[:60]), flush=True)
        f, c = parse_catch2_xml(xml_path, "focus", ctx["stems"])
        if not c and os.path.isfile(xml_path):
            findings.append(new_finding(
                "focus", "suite", "suite.focus_case_not_found", "-", case, "error",
                detail="用例名未匹配到任何用例，请核对 tools/focus_cases.txt",
                location=rel_path(args.focus_list)))
        findings += f
        cases += c
    return findings, cases


def tally(findings, key):
    """按 key 聚合成 {值: {error,warn,info,total,peers}}；peers 记录关联维度取值。"""
    result = OrderedDict()
    for item in findings:
        cell = result.setdefault(item[key], {"error": 0, "warn": 0, "info": 0,
                                             "total": 0, "peers": OrderedDict()})
        severity = item["severity"] if item["severity"] in SEVERITY_ORDER else "info"
        cell[severity] += 1
        cell["total"] += 1
        peer = item["dataset"] if key == "category" else item["category"]
        if severity != "info":
            cell["peers"][peer] = cell["peers"].get(peer, 0) + 1
    return result


def top_peers(peers, limit=3):
    ordered = sorted(peers.items(), key=lambda kv: (-kv[1], kv[0]))[:limit]
    return "、".join("%s(%d)" % (k, v) for k, v in ordered)


def finding_key(item):
    """基线对比的身份键：不含数值，避免耗时/偏差抖动被误判成“新增”。"""
    return (item["category"], item["check"], item["dataset"], item["subject"])


def apply_baseline(findings, baseline_path):
    """给当前发现打上 新增/保持 标记，并返回基线里已消除的条目。"""
    base = read_findings_csv(baseline_path)
    base_keys = set(finding_key(b) for b in base)
    cur_keys = set(finding_key(f) for f in findings)
    for item in findings:
        item["status"] = "保持" if finding_key(item) in base_keys else "新增"
    resolved = [b for b in base if finding_key(b) not in cur_keys]
    for item in resolved:
        item["status"] = "已消除"
    return resolved


def sort_findings(findings):
    order = {cid: i for i, cid in enumerate(CATEGORY_LABELS)}
    return sorted(findings, key=lambda f: (
        SEVERITY_ORDER.get(f["severity"], 3), order.get(f["category"], 99),
        f["dataset"], f["check"], f["subject"]))


DETAIL_HEADERS = ["来源", "检查类别", "检查项", "数据集/文件", "对象", "严重度",
                  "指标", "实测值", "阈值", "说明", "位置"]
DETAIL_WIDTHS = [10, 12, 26, 24, 28, 9, 14, 12, 10, 62, 34]


def detail_row(item, with_status):
    row = [SOURCE_LABELS.get(item["source"], item["source"]),
           label_of(item["category"]), item["check"], item["dataset"],
           item["subject"], item["severity"], item["metric"],
           to_number(item["value"]), to_number(item["threshold"]),
           item["detail"], item["location"]]
    if with_status:
        row.append(item.get("status", ""))
    return row


def to_number(text):
    """CSV 里的数值列转成 Excel 数字，空/非数值原样返回字符串。"""
    if text is None or text == "":
        return ""
    try:
        value = float(text)
    except (TypeError, ValueError):
        return text
    return int(value) if value == int(value) and abs(value) < 1e15 else value


def ordered_categories(by_category):
    """先按预设顺序，再追加未登记的新类别，保证 sheet 顺序稳定。"""
    known = [c for c in CATEGORY_LABELS if c in by_category]
    extra = [c for c in by_category if c not in CATEGORY_LABELS]
    return known + sorted(extra)


def add_summary_sheet(wb, findings, cases, resolved, meta, by_cat, by_ds):
    total = len(findings)
    counts = {s: sum(1 for f in findings if f["severity"] == s)
              for s in ("error", "warn", "info")}
    rows = [["报表生成时间", meta["timestamp"]],
            ["执行阶段", meta["stages"]],
            ["数据目录", meta["datas"]],
            ["构建目录", meta["build_dir"]],
            ["复用已有产物的阶段(--skip-run)", meta["skip_run"] or "无（全部实跑）"],
            ["", ""],
            ["发现总数", total],
            ["error（必须修复）", counts["error"]],
            ["warn（需确认）", counts["warn"]],
            ["info（记录/基线）", counts["info"]],
            ["涉及检查类别数", len(by_cat)],
            ["涉及文件数", len([d for d in by_ds if d != "-"])],
            ["", ""]]
    for source, label in SOURCE_LABELS.items():
        subset = [f for f in findings if f["source"] == source]
        if not subset and source not in meta["stage_set"]:
            continue
        rows.append(["%s 发现数(error/warn/info)" % label,
                     "%d (%d/%d/%d)" % (
                         len(subset),
                         sum(1 for f in subset if f["severity"] == "error"),
                         sum(1 for f in subset if f["severity"] == "warn"),
                         sum(1 for f in subset if f["severity"] == "info"))])
    if cases:
        rows += [["", ""],
                 ["用例总数", len(cases)],
                 ["用例通过", sum(1 for c in cases if c["result"] == "pass")],
                 ["用例失败", sum(1 for c in cases if c["result"] == "fail")],
                 ["用例跳过", sum(1 for c in cases if c["result"] == "skip")]]
    top_cat = sorted(by_cat.items(), key=lambda kv: -kv[1]["error"])[:5]
    top_ds = sorted(by_ds.items(), key=lambda kv: -kv[1]["error"])[:5]
    rows += [["", ""], ["error 最多的类别", top_peers(
        OrderedDict((label_of(c), v["error"]) for c, v in top_cat if v["error"]), 5)]]
    rows.append(["error 最多的文件", top_peers(
        OrderedDict((d, v["error"]) for d, v in top_ds if v["error"]), 5)])
    if meta["baseline"]:
        rows += [["", ""], ["基线文件", meta["baseline"]],
                 ["新增", sum(1 for f in findings if f.get("status") == "新增")],
                 ["保持", sum(1 for f in findings if f.get("status") == "保持")],
                 ["已消除", len(resolved)]]
    wb.add_sheet("汇总", ["项目", "值"], rows, widths=[34, 78], autofilter=False)


def _tally_style(cell):
    if cell["error"]:
        return X.STYLE_ERROR
    if cell["warn"]:
        return X.STYLE_WARN
    return X.STYLE_INFO


def add_category_sheet(wb, by_cat):
    rows, styles = [], []
    for cid in sorted(by_cat, key=lambda c: (-by_cat[c]["error"], -by_cat[c]["total"], c)):
        cell = by_cat[cid]
        rows.append([label_of(cid), cid, cell["error"], cell["warn"], cell["info"],
                     cell["total"], len([p for p in cell["peers"] if p != "-"]),
                     top_peers(cell["peers"])])
        styles.append(_tally_style(cell))
    wb.add_sheet("按类别统计",
                 ["检查类别", "类别id", "error", "warn", "info", "合计",
                  "涉及文件数", "error+warn 最多的文件(前3)"],
                 rows, widths=[16, 20, 9, 9, 9, 9, 12, 52], row_styles=styles)


def add_dataset_sheet(wb, by_ds):
    rows, styles = [], []
    for name in sorted(by_ds, key=lambda d: (-by_ds[d]["error"], -by_ds[d]["total"], d)):
        cell = by_ds[name]
        rows.append([name, cell["error"], cell["warn"], cell["info"], cell["total"],
                     len(cell["peers"]),
                     top_peers(OrderedDict((label_of(k), v)
                                           for k, v in cell["peers"].items()))])
        styles.append(_tally_style(cell))
    wb.add_sheet("按文件统计",
                 ["数据集/文件", "error", "warn", "info", "合计", "涉及类别数",
                  "error+warn 最多的类别(前3)"],
                 rows, widths=[30, 9, 9, 9, 9, 12, 52], row_styles=styles)


def add_matrix_sheet(wb, findings):
    """类别×文件矩阵，单元格是 error 条数，用于挑“先修哪一格”。"""
    cats = ordered_categories(tally(findings, "category"))
    datasets = sorted(set(f["dataset"] for f in findings))
    grid = {}
    for item in findings:
        if item["severity"] != "error":
            continue
        key = (item["dataset"], item["category"])
        grid[key] = grid.get(key, 0) + 1
    rows, styles = [], []
    for name in sorted(datasets, key=lambda d: -sum(
            grid.get((d, c), 0) for c in cats)):
        total = sum(grid.get((name, c), 0) for c in cats)
        rows.append([name] + [grid.get((name, c), 0) or "" for c in cats] + [total])
        styles.append(X.STYLE_ERROR if total else X.STYLE_INFO)
    headers = ["数据集/文件"] + [label_of(c) for c in cats] + ["合计error"]
    widths = [30] + [12] * len(cats) + [10]
    wb.add_sheet("类别x文件矩阵", headers, rows, widths=widths, row_styles=styles)


def add_cases_sheet(wb, cases):
    if not cases:
        return
    style_of = {"pass": X.STYLE_OK, "fail": X.STYLE_ERROR, "skip": X.STYLE_WARN}
    ordered = sorted(cases, key=lambda c: ({"fail": 0, "skip": 1, "pass": 2}[c["result"]],
                                           -c["duration"], c["suite"], c["case"]))
    rows = [[SOURCE_LABELS.get(c["source"], c["source"]), c["suite"], c["case"],
             c["tags"], c["result"], c["failed"], c["skips"], round(c["duration"], 3),
             label_of(c["category"]), c["dataset"], c["location"]] for c in ordered]
    wb.add_sheet("用例结果",
                 ["来源", "套件", "用例名", "标签", "结果", "失败断言数", "跳过数",
                  "耗时(s)", "检查类别", "数据集/文件", "位置"],
                 rows, widths=[10, 32, 60, 30, 8, 12, 9, 10, 14, 24, 40],
                 row_styles=[style_of[c["result"]] for c in ordered])


def add_baseline_sheet(wb, findings, resolved):
    rows, styles = [], []
    for item in resolved + [f for f in findings if f.get("status") == "新增"]:
        rows.append([item.get("status", ""), label_of(item["category"]), item["check"],
                     item["dataset"], item["subject"], item["severity"],
                     item["detail"], item["location"]])
        styles.append(X.STYLE_OK if item.get("status") == "已消除"
                      else X.style_for_severity(item["severity"]))
    if not rows:
        rows, styles = [["无变化", "", "", "", "", "", "与基线逐条一致", ""]], [X.STYLE_OK]
    wb.add_sheet("基线对比",
                 ["状态", "检查类别", "检查项", "数据集/文件", "对象", "严重度",
                  "说明", "位置"],
                 rows, widths=[10, 14, 26, 24, 28, 9, 62, 34], row_styles=styles)


def build_workbook(findings, cases, resolved, meta):
    with_status = bool(meta["baseline"])
    by_cat = tally(findings, "category")
    by_ds = tally(findings, "dataset")
    wb = X.Workbook()
    add_summary_sheet(wb, findings, cases, resolved, meta, by_cat, by_ds)
    add_category_sheet(wb, by_cat)
    add_dataset_sheet(wb, by_ds)
    add_matrix_sheet(wb, findings)

    headers = DETAIL_HEADERS + (["状态"] if with_status else [])
    widths = DETAIL_WIDTHS + ([8] if with_status else [])
    for cid in ordered_categories(by_cat):
        subset = [f for f in findings if f["category"] == cid]
        wb.add_sheet("%s(%d)" % (label_of(cid), len(subset)), headers,
                     [detail_row(f, with_status) for f in subset], widths,
                     [X.style_for_severity(f["severity"]) for f in subset])
    wb.add_sheet("全部明细", headers, [detail_row(f, with_status) for f in findings],
                 widths, [X.style_for_severity(f["severity"]) for f in findings])
    add_cases_sheet(wb, cases)
    if with_status:
        add_baseline_sheet(wb, findings, resolved)
    return wb, by_cat, by_ds


def print_tables(by_cat, by_ds):
    """终端也打一份统计，方便不开 Excel 时直接按类别/文件排修复顺序。"""
    print("\n==== 按检查类别统计（error 降序）====")
    print("%-14s %-20s %6s %6s %6s %6s  %s"
          % ("类别", "类别id", "error", "warn", "info", "合计", "主要文件"))
    for cid in sorted(by_cat, key=lambda c: (-by_cat[c]["error"], -by_cat[c]["total"], c)):
        cell = by_cat[cid]
        print("%-14s %-20s %6d %6d %6d %6d  %s"
              % (label_of(cid), cid, cell["error"], cell["warn"], cell["info"],
                 cell["total"], top_peers(cell["peers"])))
    print("\n==== 按文件统计（error 降序）====")
    print("%-32s %6s %6s %6s %6s  %s"
          % ("数据集/文件", "error", "warn", "info", "合计", "主要类别"))
    for name in sorted(by_ds, key=lambda d: (-by_ds[d]["error"], -by_ds[d]["total"], d)):
        cell = by_ds[name]
        print("%-32s %6d %6d %6d %6d  %s"
              % (name, cell["error"], cell["warn"], cell["info"], cell["total"],
                 top_peers(OrderedDict((label_of(k), v)
                                       for k, v in cell["peers"].items()))))


def parse_args(argv):
    parser = argparse.ArgumentParser(
        prog="test_report.py", description="四类测试统一 Excel 报表",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="阶段可选：audit(全量审计) regression(回归套件) focus(焦点用例) perf(性能)")
    parser.add_argument("--build-dir", default="cmake-build-release",
                        help="构建目录（需含 tests/ 下的可执行文件），默认 cmake-build-release")
    parser.add_argument("--datas", default="datas", help="数据目录，默认 datas")
    parser.add_argument("--out", default="reports", help="输出目录，默认 reports")
    parser.add_argument("--stages", default="audit,regression,focus,perf",
                        help="要执行的阶段，逗号分隔")
    parser.add_argument("--suites", default="", help="只跑指定 Catch2 套件，逗号分隔")
    parser.add_argument("--focus-list", default=os.path.join(HERE, "focus_cases.txt"),
                        help="焦点用例清单，默认 tools/focus_cases.txt")
    parser.add_argument("--repeat", type=int, default=1, help="性能阶段重复次数（取中位数）")
    parser.add_argument("--only", default="", help="审计阶段只处理文件名含该子串的数据集")
    parser.add_argument("--road-edge-clearance", type=float, default=None,
                        help="审计阶段覆盖路沿净距阈值（米）")
    parser.add_argument("--baseline", default="", help="上一轮 findings.csv，用于新增/已消除对比")
    parser.add_argument("--skip-run", nargs="?", const="all", default="",
                        metavar="阶段",
                        help="复用 out 目录已有产物而不实跑；不带值=全部阶段，"
                             "也可指定 audit,perf 等子集")
    parser.add_argument("--xlsx", default="", help="Excel 文件名，默认 测试报告_<时间戳>.xlsx")
    args = parser.parse_args(argv)
    for name in ("build_dir", "datas", "out"):
        setattr(args, name, os.path.join(ROOT, getattr(args, name))
                if not os.path.isabs(getattr(args, name)) else getattr(args, name))
    return args


def main(argv=None):
    args = parse_args(argv)
    stages = [s.strip() for s in args.stages.split(",") if s.strip()]
    unknown = [s for s in stages if s not in ("audit", "regression", "focus", "perf")]
    if unknown:
        print("未知阶段: %s" % ",".join(unknown), file=sys.stderr)
        return 2

    ctx = {"out": args.out,
           "logs": os.path.join(args.out, "logs"),
           "xml": os.path.join(args.out, "xml"),
           "stems": dataset_stems(args.datas)}
    for path in (ctx["out"], ctx["logs"], ctx["xml"]):
        if not os.path.isdir(path):
            os.makedirs(path)

    findings, cases = [], []
    print("执行阶段: %s；复用已有产物: %s"
          % (",".join(stages), args.skip_run or "无"), flush=True)
    if "audit" in stages:
        findings += stage_audit(args, ctx)
    if "regression" in stages:
        f, c = stage_regression(args, ctx)
        findings += f
        cases += c
    if "focus" in stages:
        f, c = stage_focus(args, ctx)
        findings += f
        cases += c
    if "perf" in stages:
        findings += stage_perf(args, ctx)

    findings = sort_findings(findings)
    resolved = apply_baseline(findings, args.baseline) if args.baseline else []
    meta = {"timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "stages": ",".join(stages), "stage_set": set(stages),
            "datas": rel_path(args.datas), "build_dir": rel_path(args.build_dir),
            "skip_run": args.skip_run,
            "baseline": rel_path(args.baseline) if args.baseline else ""}

    csv_path = os.path.join(args.out, "findings.csv")
    write_findings_csv(csv_path, findings)
    workbook, by_cat, by_ds = build_workbook(findings, cases, resolved, meta)
    xlsx_name = args.xlsx or time.strftime("测试报告_%Y%m%d_%H%M%S.xlsx")
    xlsx_path = os.path.join(args.out, xlsx_name)
    workbook.save(xlsx_path)

    print_tables(by_cat, by_ds)
    errors = sum(1 for f in findings if f["severity"] == "error")
    warns = sum(1 for f in findings if f["severity"] == "warn")
    print("\n发现 %d 条（error %d / warn %d / info %d），用例 %d 条"
          % (len(findings), errors, warns, len(findings) - errors - warns, len(cases)))
    print("Excel : %s" % rel_path(xlsx_path))
    print("CSV   : %s" % rel_path(csv_path))
    # 退出码只反映“报表本身是否可信”：工具级错误（可执行文件缺失、异常退出、XML
    # 解析失败、焦点用例名对不上）意味着数据不全；用例断言失败属于被测对象的问题，
    # 已经体现在统计表里，不应让报表命令失败，否则 CI 无法区分二者。
    tool_errors = [f for f in findings
                   if f["check"].startswith("suite.") and f["severity"] == "error"]
    if tool_errors:
        print("\n注意：%d 条工具级错误，报表可能不完整：" % len(tool_errors),
              file=sys.stderr)
        for f in tool_errors:
            print("  %s %s %s" % (f["check"], f["subject"], f["detail"]),
                  file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
