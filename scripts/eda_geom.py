# -*- coding: utf-8 -*-
"""EasyEDA 原生几何 source 的解析（正确版）。

source 是 EasyEDA 原生的「complex polygon」编码，语法：

    ['R', x0, y0, w, h, rot, 0]        矩形（只在禁布区里出现）
    [x0, y0, 'L', x1, y1, ..., 'ARC', sweep, ex, ey, 'L', ...]
                                       折线 + 圆弧；外层可能是多个轮廓的数组

关键事实（逐条用真实数据验证过，见 validate.py）：
  * ARC 后面跟 3 个数：**sweep 角度、终点 x、终点 y**。
    sweep 是角度（可为负），(ex,ey) 是圆弧终点 —— 注意终点也是路径顶点。
    ⚠ 早期版本把后两个数当圆心，导致圆弧被当成孤立点丢掉或算错。
  * 一个轮廓若只有两个 ARC（如 `P0,'ARC',180,e1,'ARC',180,P0`），
    这是两个半圆拼成的整圆。180° 圆弧的圆心就是弦中点。
  * 非 180° 圆弧的圆心由「起点、终点、扫角」决定：弦长 L 与半径 r 满足
    L = 2 r sin(|θ|/2)；圆心在弦中垂线上，sweep>0 表示逆时针。
  * 外层数组的每个元素是一个**独立轮廓**。铺铜用 even-odd 填充规则：
    点被奇数个轮廓包住才算在铜里（这样内部挖空/禁布孔才正确）。
"""
import math


def is_num(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _arc_points(start, end, sweep_deg, steps):
    """把一个圆弧起点->终点按 sweep 展开成中间点（不含起始点，含终点前）。"""
    sx, sy = start
    ex, ey = end
    if abs(sweep_deg) < 1e-9:
        return []
    chord = math.hypot(ex - sx, ey - sy)
    if chord < 1e-9:
        return []
    half = math.radians(abs(sweep_deg)) / 2.0
    sin_half = math.sin(half)
    if sin_half < 1e-12:
        return []
    r = (chord / 2.0) / sin_half
    mx, my = (sx + ex) / 2.0, (sy + ey) / 2.0
    if abs(abs(sweep_deg) - 180.0) < 1e-6:
        cx, cy = mx, my
    else:
        h = math.sqrt(max(r * r - (chord / 2.0) ** 2, 0.0))
        ux, uy = (ex - sx) / chord, (ey - sy) / chord
        # 逆时针时圆心在弦方向左侧
        if sweep_deg > 0:
            nx, ny = -uy, ux
        else:
            nx, ny = uy, -ux
        cx, cy = mx + h * nx, my + h * ny
    a0 = math.atan2(sy - cy, sx - cx)
    a1 = math.atan2(ey - cy, ex - cx)
    if sweep_deg > 0:
        da = (a1 - a0) % (2 * math.pi)
    else:
        da = -((a0 - a1) % (2 * math.pi))
    out = []
    for k in range(1, steps):
        t = k / float(steps)
        ang = a0 + da * t
        out.append((cx + r * math.cos(ang), cy + r * math.sin(ang)))
    return out


def parse_contour(seg, arc_steps=32):
    """解析单个轮廓成顶点列表。"""
    pts = []
    i = 0
    n = len(seg)
    while i < n:
        v = seg[i]
        if isinstance(v, str):
            cmd = v.upper()
            if cmd in ("ARC", "CARC"):
                if i + 3 >= n:
                    break
                sweep = float(seg[i + 1])
                ex, ey = float(seg[i + 2]), float(seg[i + 3])
                if pts:
                    pts.extend(_arc_points(pts[-1], (ex, ey), sweep, arc_steps))
                pts.append((ex, ey))
                i += 4
                continue
            # L / 其它命令标号，跳过
            i += 1
            continue
        if not is_num(v):
            i += 1
            continue
        if i + 1 < n and is_num(seg[i + 1]):
            pts.append((float(v), float(seg[i + 1])))
            i += 2
        else:
            i += 1
    return pts


def parse_rect(source, arc_steps=32):
    """['R', x0, y0, w, h, rot, 0] -> 四角。

    局部角点 (0,0),(w,0),(w,-h),(0,-h)，绕 (x0,y0) 逆时针转 rot。
    已验证：rot=0→x[x0,x0+w] y[y0-h,y0]；rot=90→x[x0,x0+h] y[y0,y0+w]；
            rot=180→x[x0-w,x0] y[y0,y0+h]；rot=270→x[x0-h,x0] y[y0-w,y0]。
    """
    f = [float(v) for v in source[1:6] if is_num(v)]
    x0, y0, w, h = f[0], f[1], f[2], f[3]
    rot = f[4] if len(f) > 4 else 0.0
    ca, sa = math.cos(math.radians(rot)), math.sin(math.radians(rot))
    return [(x0 + lx * ca - ly * sa, y0 + lx * sa + ly * ca)
            for lx, ly in ((0.0, 0.0), (w, 0.0), (w, -h), (0.0, -h))]


def parse_contours(source, arc_steps=32):
    """解析成轮廓列表。

    返回 list[list[(x,y)]]。矩形返回单个四角轮廓；多轮廓数组逐个解析。
    """
    if not source:
        return []
    if isinstance(source, list) and source and source[0] == "R":
        return [parse_rect(source, arc_steps)]
    segs = source if isinstance(source[0], list) else [source]
    out = []
    for seg in segs:
        pts = parse_contour(seg, arc_steps)
        if len(pts) >= 3:
            out.append(pts)
    return out


def bbox_of(contours):
    xs = [p[0] for c in contours for p in c]
    ys = [p[1] for c in contours for p in c]
    if not xs:
        return None
    return {"minX": min(xs), "maxX": max(xs), "minY": min(ys), "maxY": max(ys)}


def point_in_contour(pt, poly):
    if len(poly) < 3:
        return False
    x, y = pt
    inside = False
    n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        if (y1 > y) != (y2 > y):
            xin = (x2 - x1) * (y - y1) / (y2 - y1) + x1
            if x < xin:
                inside = not inside
    return inside


def point_in_complex(pt, contours):
    """even-odd：被奇数个轮廓包住才算在内（内部挖空会被正确排除）。"""
    hits = 0
    for c in contours:
        if point_in_contour(pt, c):
            hits += 1
    return hits % 2 == 1
