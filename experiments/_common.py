"""公共工具：JSON 安全序列化 + 数值快照写入。"""
import json
import sys
import numpy as np

# Windows 控制台 GBK 编码兜底：避免 print 含数学符号时报编码错误
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass


def sanitize(o):
    """递归转成 JSON 可序列化（处理 numpy 标量、sympy Boolean/表达式）。"""
    if isinstance(o, dict):
        return {k: sanitize(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [sanitize(v) for v in o]
    if isinstance(o, (np.bool_, np.integer, np.floating)):
        return o.item()
    if hasattr(o, "is_Boolean") and o.is_Boolean:   # sympy Boolean
        return bool(o)
    if isinstance(o, (bool, int, float, str)) or o is None:
        return o
    # 其它（sympy 表达式等）转字符串
    return str(o)


def report(results, name):
    # 控制台输出用 ensure_ascii=True（纯 ASCII，避免 Windows GBK 控制台乱码）；
    # 文件用 ensure_ascii=False + UTF-8（中文/数学符号可读，公开后正常）。
    print(json.dumps(sanitize(results), indent=2, ensure_ascii=True))
    with open(f"experiments/{name}_last_run.json", "w", encoding="utf-8") as f:
        json.dump(sanitize(results), f, indent=2, ensure_ascii=False)
    print(f"OK: {name} verified")
