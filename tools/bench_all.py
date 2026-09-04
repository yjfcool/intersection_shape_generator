#!/usr/bin/env python3
"""逐份数据集测量单路口生成耗时，带单例超时保护。

diag_gen_time 一次性跑完全部数据时，任何一份病态数据都会拖垮整轮测量，
且 stdout 被重定向后是全缓冲的，进程被杀时看不到任何已完成的结果。
本脚本改为「一份数据一个进程」，逐份落盘，超时单独标记 TIMEOUT。

用法:
  tools/bench_all.py [--bin 可执行文件] [--dir 数据目录] [--timeout 秒]
                     [--repeat 次数] [--out 结果.csv] [--only 子串,...]
"""
import argparse
import csv
import glob
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HARD_LIMIT_MS = 15000.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=os.path.join(
        ROOT, "cmake-build-release", "tests", "diag_gen_time"))
    ap.add_argument("--dir", default=os.path.join(ROOT, "datas"))
    ap.add_argument("--timeout", type=float, default=120.0)
    ap.add_argument("--repeat", type=int, default=1)
    ap.add_argument("--out", default=os.path.join(ROOT, "bench_result.csv"))
    ap.add_argument("--only", default="")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.dir, "*.json")))
    if args.only:
        keys = [k for k in args.only.split(",") if k]
        files = [f for f in files if any(k in os.path.basename(f) for k in keys)]
    if not files:
        print("no dataset found", file=sys.stderr)
        return 2

    rows = []
    for path in files:
        name = os.path.basename(path)
        cmd = [args.bin, path, "--repeat", str(args.repeat)]
        t0 = time.time()
        status, ms = "ok", None
        try:
            p = subprocess.run(cmd, capture_output=True, text=True,
                               timeout=args.timeout)
            wall = (time.time() - t0) * 1000.0
            out = p.stdout.strip().splitlines()
            for line in out:
                if line.startswith(path) and line.endswith("ms"):
                    parts = line.split()
                    ms = float(parts[-2])
                    if parts[-3] == "ok=0":
                        status = "gen_false"
            if ms is None:
                status, ms = "crash", wall
        except subprocess.TimeoutExpired:
            status, ms = "TIMEOUT", args.timeout * 1000.0

        flag = "  <-- OVER 15s" if ms > HARD_LIMIT_MS or status == "TIMEOUT" else ""
        print("%-28s %-9s %10.1f ms%s" % (name, status, ms, flag), flush=True)
        rows.append({"dataset": name, "status": status, "ms": round(ms, 1)})

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["dataset", "status", "ms"])
        w.writeheader()
        w.writerows(rows)

    over = [r for r in rows if r["ms"] > HARD_LIMIT_MS or r["status"] == "TIMEOUT"]
    total = sum(r["ms"] for r in rows)
    print("\n---- %d datasets, total %.1f s, over-15s: %d ----"
          % (len(rows), total / 1000.0, len(over)))
    for r in sorted(over, key=lambda r: -r["ms"]):
        print("  %-28s %10.1f ms (%s)" % (r["dataset"], r["ms"], r["status"]))
    print("result -> %s" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
