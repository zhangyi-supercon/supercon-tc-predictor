# -*- coding: utf-8 -*-
"""
supercon_formula_lib.py — 自包含 SuperCon 配方特征工程库 v1.1
替代缺失的 R205_supercon_formula_predictor / R205_元素属性表 依赖。
功能:
  1. parse_formula(formula) -> [(元素, 化学计量), ...]
  2. build_features(formula) -> dict 81 特征 (Hamidieh 2018 特征工程, 与 train.csv 列序一致)
  3. ELEMENT_PROPS: 元素属性表 (逆特征工程反推, 与训练数据 MAE=0)
"""
import re
import math

# ===== 元素属性表 =====
# 2026-08-16 逆特征工程: 从 train.csv(81特征) + unique_m.csv(元素计数) 最小二乘反推,
# 与训练特征完全一致 (验证 MAE=0)。生成脚本: tmp_infer_props_v94.py
from supercon_element_props_data import ELEMENT_PROPS

# ===== 化学式解析 =====
# 支持: MgB2, YBa2Cu3O7, Ba0.2La1.8Cu1O4, FeTe0.5Se0.5, Bi-2212, (LaH10)2, LaH10
_ELEM_RE = re.compile(r'([A-Z][a-z]?)([0-9]*\.?[0-9]*)')


def parse_formula(formula):
    """解析化学式 -> [(元素符号, 化学计量), ...]。支持小数系数/括号组/常见别名。"""
    f = formula.strip()
    if not f:
        raise ValueError('空化学式')
    # 常见别名
    alias = {
        'bi-2212': 'Bi2Sr2CaCu2O8', 'bi2212': 'Bi2Sr2CaCu2O8',
        'bi-2223': 'Bi2Sr2Ca2Cu3O10', 'bi2223': 'Bi2Sr2Ca2Cu3O10',
        'ybco': 'YBa2Cu3O7', 'bscco': 'Bi2Sr2CaCu2O8',
        'nb3sn': 'Nb3Sn', 'mgb2': 'MgB2', 'fese': 'FeSe',
        'lah10': 'LaH10', 'h3s': 'H3S', 'lh10': 'LaH10',
    }
    key = f.lower().replace(' ', '')
    if key in alias:
        f = alias[key]

    # 展开括号组: (A2B3)4 -> A8B12
    while '(' in f:
        m = re.search(r'\(([^()]*)\)([0-9]*\.?[0-9]*)', f)
        if not m:
            raise ValueError('括号不匹配: %s' % formula)
        inner, mult = m.group(1), float(m.group(2) or 1)
        inner_comp = _parse_plain(inner)
        expanded = ''.join('%s%s' % (el, cnt * mult) for el, cnt in inner_comp)
        f = f[:m.start()] + expanded + f[m.end():]

    comp = _parse_plain(f)
    if not comp:
        raise ValueError('无法解析化学式: %s' % formula)
    # 合并同元素
    merged = {}
    for el, cnt in comp:
        merged[el] = merged.get(el, 0.0) + cnt
    return [(el, cnt) for el, cnt in merged.items() if cnt > 0]


def _parse_plain(part):
    comp = []
    for m in _ELEM_RE.finditer(part):
        el = m.group(1)
        cnt = float(m.group(2)) if m.group(2) else 1.0
        comp.append((el, cnt))
    return comp


# ===== 81 维特征工程 (Hamidieh 2018) =====
_PROP_KEYS = ['atomic_mass', 'fie', 'atomic_radius', 'Density',
              'ElectronAffinity', 'FusionHeat', 'ThermalConductivity', 'Valence']
_STATS = ['mean', 'wtd_mean', 'gmean', 'wtd_gmean', 'entropy', 'wtd_entropy',
          'range', 'wtd_range', 'std', 'wtd_std']

FEATURE_COLS = ['number_of_elements'] + ['%s_%s' % (s, p) for p in _PROP_KEYS for s in _STATS]


def _stats(prop_vals, weights):
    """对一组属性值计算10种统计量 (Hamidieh 2018 / magpie 定义, 已用真实数据反推校准)。
    定义: mean=简单平均; wtd_mean=按化学计量加权平均; gmean/wtd_gmean=几何平均;
    entropy: p_i=v_i/Σv; wtd_entropy: p_i=w_i*v_i/Σ(w*v); range=max-min;
    wtd_range=max(v*p)-min(v*p) p=w/Σw; std/wtd_std=标准差。"""
    n = len(prop_vals)
    if n == 0:
        return [0.0] * 10
    w = [max(x, 1e-9) for x in weights]
    wsum = sum(w)
    mean = sum(prop_vals) / n
    wtd_mean = sum(v * wi for v, wi in zip(prop_vals, w)) / wsum
    pos = [max(v, 1e-9) for v in prop_vals]
    gmean = math.exp(sum(math.log(v) for v in pos) / n)
    wgmean = math.exp(sum(math.log(v) * wi for v, wi in zip(pos, w)) / wsum)
    p = [x / wsum for x in w]
    # entropy: 属性值归一化比例
    vsum = sum(prop_vals)
    pv = [x / vsum for x in prop_vals] if vsum != 0 else [0.0] * n
    entropy = -sum(x * math.log(x) for x in pv if x > 0)
    # wtd_entropy: 加权属性值归一化
    wv = [wi * vi for wi, vi in zip(w, prop_vals)]
    wvsum = sum(wv)
    pwv = [x / wvsum for x in wv] if wvsum != 0 else [0.0] * n
    wtd_entropy = -sum(x * math.log(x) for x in pwv if x > 0)
    rng = max(prop_vals) - min(prop_vals)
    vp = [vi * pi for vi, pi in zip(prop_vals, p)]
    wtd_rng = max(vp) - min(vp)
    std = math.sqrt(sum((v - mean) ** 2 for v in prop_vals) / n) if n > 1 else 0.0
    wtd_std = math.sqrt(sum(w[i] * (prop_vals[i] - wtd_mean) ** 2 for i in range(n)) / wsum) if n > 1 else 0.0
    return [mean, wtd_mean, gmean, wgmean, entropy, wtd_entropy, rng, wtd_rng, std, wtd_std]


def build_features(formula):
    """化学式 -> 81维特征 dict (与 train.csv 列序一致)。未知元素抛 ValueError。"""
    comp = parse_formula(formula)
    weights = [c for _, c in comp]
    feats = {'number_of_elements': float(len(comp))}
    for prop in _PROP_KEYS:
        vals = []
        for el, _ in comp:
            if el not in ELEMENT_PROPS:
                raise ValueError('未知元素: %s (%s)' % (el, formula))
            idx = _PROP_KEYS.index(prop)
            vals.append(ELEMENT_PROPS[el][idx])
        for s, v in zip(_STATS, _stats(vals, weights)):
            feats['%s_%s' % (s, prop)] = v
    return feats


if __name__ == '__main__':
    import sys, io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    for f in ['MgB2', 'YBa2Cu3O7', 'Nb3Sn', 'FeSe', 'Ba0.2La1.8Cu1O4', 'LaH10', 'H3S', 'Bi-2212']:
        try:
            comp = parse_formula(f)
            feats = build_features(f)
            print('%-16s %s -> %d feats' % (f, comp, len(feats)))
        except Exception as e:
            print('%-16s ERROR: %s' % (f, e))
