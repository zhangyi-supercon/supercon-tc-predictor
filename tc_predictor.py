# -*- coding: utf-8 -*-
"""超导Tc预测工具 - 核心: 化学式解析 + 预测 + 边界标注(R138 B1-B6) + CLI
用法:
  python tc_predictor.py "YBa2Cu3O7"
  python tc_predictor.py "Hg0.7Pb0.3Ba2Ca2Cu3O8" --detail
支持: 元素符号(大小写混合) + 原子数(整数/小数), 如 Ba0.2La1.8Cu1O4, H2S1
边界标注: 读取 boundary_rules.json (源自 R138 外推校准, 规则 B1-B6),
          每次预测输出信任区(green/yellow/red) + 警告列表, 不改变预测数值。
"""
import os, sys, io, re, json
import numpy as np
try:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
except Exception:
    pass

TOOL_DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_PATH = os.path.join(TOOL_DIR, 'tc_model.joblib')
BOUNDARY_PATH = os.path.join(TOOL_DIR, 'boundary_rules.json')

# 元素周期表(按unique_m列顺序)
ELEMENTS = ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn']
ELEM_SET = set(ELEMENTS)

# ---------- 81维特征桥接 (2026-08-16 第96帧: 反推表反哺常压轨) ----------
# 第94帧逆特征工程产出 supercon_formula_lib (ELEMENT_PROPS MAE=0),
# 使部署侧能真实计算 Hamidieh 81 维特征 -> boundary_rules design_rules (r3/B7)
# 从"含B就提示"升级为"真实特征区间判定"。桥接失败自动降级为元素级提示。

# 常用材料别名 -> 标准化学式 (2026-08-16 第97帧: 别名误解析缺陷修复)
# 背景: YBCO 被解析为 Y-B-C-O(含硼) -> 预测16.3K(真实85.5K,偏差-69K) + r3误报;
# BSCCO -> B-S-C-C-O -> 6.3K(真实78K)。工程材料通用缩写必须归一化。
FORMULA_ALIASES = {
    'YBCO': 'YBa2Cu3O7',
    'Y123': 'YBa2Cu3O7',
    'Y-123': 'YBa2Cu3O7',
    'BSCCO': 'Bi2Sr2CaCu2O8',
    'BI2212': 'Bi2Sr2CaCu2O8',
    'BI-2212': 'Bi2Sr2CaCu2O8',
    'BI2223': 'Bi2Sr2Ca2Cu3O10',
    'BI-2223': 'Bi2Sr2Ca2Cu3O10',
    'NB3SN': 'Nb3Sn',
    'MGB2': 'MgB2',
    'FESE': 'FeSe',
    'NBAL': 'Nb3Al',
    'YB6': 'YB6',
}

def normalize_formula(formula):
    """别名/缩写 -> 标准化学式。大小写不敏感; 无匹配原样返回。"""
    if not formula:
        return formula
    key = formula.strip().upper().replace(' ', '')
    return FORMULA_ALIASES.get(key, formula.strip())

_supercon_lib = None
_supercon_lib_attempted = False


def _get_supercon_lib():
    """懒加载 supercon_formula_lib (第94帧反推特征库)。失败返回 None 不抛异常。"""
    global _supercon_lib, _supercon_lib_attempted
    if _supercon_lib_attempted:
        return _supercon_lib
    _supercon_lib_attempted = True
    try:
        sw_dir = os.path.join(TOOL_DIR, '..', '工具', 'supercon_web')
        if sw_dir not in sys.path:
            sys.path.insert(0, sw_dir)
        import supercon_formula_lib
        _supercon_lib = supercon_formula_lib
    except Exception:
        _supercon_lib = None
    return _supercon_lib


def _build_design_feats(formula):
    """化学式 -> 81维特征 dict; 失败返回 None (桥接缺失/未知元素/解析失败)。"""
    lib = _get_supercon_lib()
    if lib is None:
        return None
    try:
        return lib.build_features(formula)
    except Exception:
        return None


def _parse_interval(s):
    """(155, inf] -> (lo, hi, lo_open, hi_open)。解析失败返回 None。"""
    m = re.match(r'^[\[\(]([^,]+),\s*([^\]\)]+)[\]\)]$', s.strip())
    if not m:
        return None
    lo_s, hi_s = m.group(1).strip(), m.group(2).strip()

    def num(x):
        if x in ('inf', '∞'):
            return float('inf')
        if x in ('-inf', '-∞'):
            return float('-inf')
        return float(x)

    return (num(lo_s), num(hi_s), s.strip().startswith('('), s.strip().endswith(')'))


def _in_interval(v, iv):
    if iv is None:
        return False
    lo, hi, lo_open, hi_open = iv
    if lo_open:
        if not (v > lo):
            return False
    else:
        if not (v >= lo):
            return False
    if hi_open:
        if not (v < hi):
            return False
    else:
        if not (v <= hi):
            return False
    return True


def _eval_design_rule(rule, counts, feats):
    """判定一条 design_rules 规则是否命中: 所有 feature 区间同时满足。
    feature 名=元素符号 -> counts 计数; 否则 -> 81维特征值。
    无 features 字段的规则(如 composite pattern 结构)不参与特征判定=恒不命中,
    它们由各自独立逻辑处理(R152-composite 等)。"""
    feats_list = rule.get('features')
    if not feats_list:
        return False
    for fspec in feats_list:
        fname = fspec.get('feature')
        iv = _parse_interval(fspec.get('range', ''))
        if iv is None:
            return False
        if fname in ELEM_SET:
            v = counts.get(fname, 0.0)
        elif feats is not None and fname in feats:
            v = feats[fname]
        else:
            return False
        if not _in_interval(v, iv):
            return False
    return True

# ---------- R138 边界标注 (B1-B6) ----------
_boundary_cache = None

def _load_boundary():
    """懒加载 boundary_rules.json; 缺失时返回空 dict(标注降级为仅词表检查)"""
    global _boundary_cache
    if _boundary_cache is None:
        try:
            with open(BOUNDARY_PATH, encoding='utf-8') as f:
                _boundary_cache = json.load(f)
        except Exception:
            _boundary_cache = {}
    return _boundary_cache

# R161 (2026-08-19 第161帧): B10 辅助——other族候选是否含高Debye温度元素载体
# 证据: Choudhary & Garrity npj Comput Mater 2022 Fig.1c——高thetaD(>300K)概率高的元素
#       = 1s/2p/3p轻元素(强共价键)与轻3d过渡金属; 重过渡金属与弱键元素 thetaD 低
HIGH_THETAD_ELEMENTS = frozenset([
    # 1s/2p/3p 轻元素 (强共价键, thetaD>300K 概率高)
    'B','C','N','O','F','Mg','Al','Si','P','S','Cl',
    # 轻3d过渡金属 (文献 Fig.1c: Sc/Ti/V/Cr/Mn/Fe/Co/Ni/Cu/Zn 高thetaD概率显著)
    'Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
])


def _has_high_thetaD_elements(counts):
    """other族候选是否含高Debye温度元素载体(BCS高Tc机理前置)。
    空counts返回False(保守: 无法判定即视为缺失)。"""
    if not counts:
        return False
    return any(e in HIGH_THETAD_ELEMENTS for e in counts)


def detect_family(counts):
    """启发式家族判定: hydride(氢化物) / cuprate(铜氧化物) / iron(铁基) / other
    R152: 升级为复合标签——骨架+间隔层+掺杂 (R151 任务3结论: 复合标签解释方差69.8%)
    如 cuprate_Hg_Ag = 铜氧骨架 + Hg间隔层 + Ag掺杂
    """
    if not counts:
        return 'other'
    total = sum(counts.values())
    h_frac = counts.get('H', 0.0) / total
    if h_frac >= 0.5:
        return 'hydride'
    fam = 'other'
    if 'Cu' in counts and 'O' in counts:
        fam = 'cuprate'
    elif 'Fe' in counts and any(s in counts for s in ('As', 'Se', 'P')):
        fam = 'iron'
    if fam == 'other':
        return fam
    # 复合标签: 间隔层元素(骨架阳离子) + 掺杂元素(非骨架主元)
    spacer = None
    for s in ('Hg', 'Tl', 'Bi', 'Pb', 'Ba', 'Sr'):
        if counts.get(s, 0) > 0:
            spacer = s
            break
    dopants = [e for e in ('Ag', 'Mg', 'Cd', 'Al', 'Ga', 'B', 'Zn', 'Ni') if counts.get(e, 0) > 0]
    tag = fam
    if spacer:
        tag += '_' + spacer
    if dopants:
        tag += '_' + '+'.join(dopants)
    return tag

def annotate_boundary(formula, counts, tc_pred):
    """按 R138 规则 B1-B6 标注预测输出。返回 dict; 不改动 tc_pred 数值。
    规则源: tc_predictor/boundary_rules.json (R139 从 r138_designer_boundary.json 提取)
    """
    b = _load_boundary()
    vocab = b.get('model_vocab') or ELEMENTS
    rules = b.get('rules') or {}
    thr = b.get('zone_thresholds') or {'green_max_tc': 80.0, 'yellow_min_tc': 80.0,
                                       'yellow_max_tc': 100.0, 'red_min_tc': 100.0}
    red_min = thr.get('red_min_tc', 100.0)
    yellow_min = thr.get('yellow_min_tc', 80.0)

    fam = detect_family(counts)
    unknown = [e for e in counts if e not in vocab]

    # 信任区: red = 高Tc 或 氢化物; yellow = 80~100K; green = 其余
    if tc_pred >= red_min or fam == 'hydride':
        zone = 'red'
    elif tc_pred >= yellow_min:
        zone = 'yellow'
    else:
        zone = 'green'

    warnings = []
    if tc_pred >= red_min:
        warnings.append({'id': 'B1', 'text': rules.get('B1', '预测>=100K: 外推区, 需文献验证, 不得作为设计结论')})
        warnings.append({'id': 'B4', 'text': rules.get('B4', '禁止在真实Tc>=100K区域用模型排序(外推区反相关)')})
    if fam == 'hydride':
        warnings.append({'id': 'B3', 'text': rules.get('B3', '氢化物族: 模型预测仅为下界参考, 真实Tc可能远超输出')})
    # R156 (2026-08-19 第156帧): B8——外部文献注入的高Tc低估风险标注
    # 证据: Xie et al. npj Comput Mater 2022 (Eliashberg/Allen-Dynes公式对高Tc常规超导体系统性低估)
    # 非氢化物 + 预测>=80K(yellow区起点) -> 附加低估风险提示, 预测视为下界参考
    elif tc_pred >= yellow_min and fam != 'other':
        warnings.append({'id': 'B8', 'text': rules.get('B8', '高Tc常规超导候选(非氢化物, 预测>=80K): 文献证据表明统计模型在高Tc端可能低估——预测应视为下界参考, 高Tc设计结论需第一性原理/实验复核')})
    # R159 (2026-08-19 第159帧): B9——正向物理锚点: 常规BCS超导的Tc物理上限
    # 证据: Choudhary & Garrity npj Comput Mater 2022——BCS预筛选(高thetaD+高N0)+DFT-PT+McMillan-Allen-Dynes,
    #       1058材料筛出105个Tc>=5K常规超导, 代表材料MoN/VC/VTe/KB6等; 常压常规超导实证上限MgB2=39K
    # 语义: other族(非铜氧/非铁基/非氢化物)=常规BCS类, 预测>=40K 与BCS物理上限冲突=伪相关风险
    # 与B8互补: B8管"高Tc被低估"(下界参考), B9管"other族高Tc可能根本不存在"(伪相关否定)
    elif fam == 'other' and tc_pred >= 40.0:
        warnings.append({'id': 'B9', 'text': rules.get('B9', '非铜氧/非铁基/非氢化物(常规BCS类)预测>=40K: 常压常规超导实证上限为MgB2=39K (Choudhary & Garrity npj Comput Mater 2022)——该预测与BCS物理上限冲突, 判定伪相关, 仅作排序线索不作结论, 需实验/文献复核')})

    # R161 (2026-08-19 第161帧): B10——机理前置否定: iBCS预筛第一关缺失载体
    # 证据: Choudhary & Garrity npj Comput Mater 2022——iBCS预筛=高Debye温度(thetaD>300K)+高DOS(N(0)>1),
    #       17419材料筛出1736; 高thetaD元素=1s/2p/3p轻元素与轻3d过渡金属(强共价键)
    # 语义: other族高Tc预测先过B9(结果层否定), 再过B10(机理前置层否定)——合起来=BCS双路监督。
    #       无高thetaD元素体系即使预测20~39K(小于B9阈值), 也没有常规BCS高Tc物理载体
    # 独立 if (非elif): B9(结果层) 与 B10(机理层) 可同时触发, 双路监督不互斥
    if fam == 'other' and tc_pred >= 30.0 and not _has_high_thetaD_elements(counts):
        warnings.append({'id': 'B10', 'text': rules.get('B10', 'other族(常规BCS类)候选不含高德拜温度元素: iBCS预筛第一关(Choudhary & Garrity npj Comput Mater 2022)要求高Debye温度(thetaD>300K)与高N(0)>1——缺机理前置条件=BCS高Tc不可能, 预测>=30K高度可疑(伪相关), 仅作排序线索不作结论, 需实验/文献复核')})

    if unknown:
        warnings.append({'id': 'B5', 'text': rules.get('B5', '超出模型词表: %s, 模型无法表示该体系' % ', '.join(unknown))})
    if zone != 'green':
        warnings.append({'id': 'B6', 'text': rules.get('B6', '高Tc候选必须附压力条件字段, 缺压力=不完整设计(压力是氢化物Tc第一变量)')})

    # R152: r3' 规则接入部署链路 (boundary_rules v3 design_rules.r3)
    # B系高Tc岛: 含B AND range_ThermalConductivity>155 AND wtd_range_ThermalConductivity>66.3
    # 部署侧仅有化学式(无法算热导跨度特征), 故按元素级信号标注: 含B → 提示 r3 候选
    design = b.get('design_rules') or {}
    # R156-v2 (2026-08-16 第96帧): 81维特征桥接——design_rules 从"元素级提示"升级为"真实特征区间判定"
    feats = _build_design_feats(formula)
    if feats is not None:
        hit_rules = []
        for rid, rule in design.items():
            if not isinstance(rule, dict):
                continue
            if rule.get('type') not in ('island_positive', 'island_positive_sparse'):
                continue
            if _eval_design_rule(rule, counts, feats):
                hit_rules.append(rid)
        if hit_rules:
            warnings.append({'id': 'design', 'text': '特征规则命中(%s): %s (81维特征真实判定, 边界规则v3)' % (
                ','.join(hit_rules), '; '.join(str(design.get(r, {}).get('name', r)) for r in hit_rules))})
    else:
        # 降级路径: 桥接不可用 -> 保留旧元素级提示 (含B)
        r3 = design.get('r3')
        if r3 and counts.get('B', 0) > 0:
            warnings.append({'id': 'r3', 'text': 'r3规则命中(含B, 特征库不可用降级): B系高Tc岛候选——完整判定需 range_ThermalConductivity>155 且 wtd_range_ThermalConductivity>66.3 (边界规则v3), 请用数据集特征核验'})
    # 复合标签识别 (R152 任务3): 规则层可识别 Hg-La-Ag 类骨架+间隔层+掺杂
    if fam != 'hydride' and fam != 'other':
        warnings.append({'id': 'R152-composite', 'text': '复合标签: %s (骨架_间隔层_掺杂, R151结论: 复合标签解释方差69.8%%)' % fam})

    return {
        'trust_zone': zone,
        'family': fam,
        'warnings': warnings,
        'pressure_boundary': '常压模型: 预测仅对常压体系可信; 加压/氢化物体系需实验压力条件验证 (R138 B2)',
        'vocab_unknown': unknown,
    }

# ---------- 化学式解析 ----------
def parse_formula(formula):
    """解析化学式 -> {元素: 原子数}。支持 Ba0.2La1.8Cu1O4 / H2S1 / YBa2Cu3O7
    规则: 大写字母开头+可选小写=元素符号, 后跟可选数字(整数/小数), 无数字=1
    返回 (counts, error)
    """
    f = formula.strip().replace(' ', '')
    if not f:
        return None, '空化学式'
    tokens = re.findall(r'([A-Z][a-z]?)(\d*\.?\d*)?', f)
    joined = ''.join(t + n for t, n in tokens)
    if joined != f:
        return None, '无法解析(包含无法识别的字符): %s' % f
    counts = {}
    for sym, num in tokens:
        if sym not in ELEM_SET:
            return None, '未知元素: %s' % sym
        n = float(num) if num else 1.0
        counts[sym] = counts.get(sym, 0.0) + n
    return counts, None

def counts_to_vector(counts):
    """元素计数 -> 86维比例向量(按ELEMENTS顺序, 归一化)"""
    v = np.zeros(len(ELEMENTS))
    for sym, n in counts.items():
        v[ELEMENTS.index(sym)] = n
    s = v.sum()
    return v / s if s > 0 else v

class TcPredictor:
    def __init__(self, model_path=MODEL_PATH):
        d = joblib_load(model_path)
        self.model = d['model']
        self.meta = d['meta']
        self.elements = self.meta['elements']

    def predict_formula(self, formula):
        counts, err = parse_formula(formula)
        if err:
            return None, err
        vec = counts_to_vector(counts)
        tc = float(self.model.predict([vec])[0])
        ann = annotate_boundary(formula, counts, tc)
        warnings = list(ann['warnings'])
        # v0.1.1 (2026-08-26 文献锚点探索): 高Tc金属氮化物外推区提示
        # 证据: 训练集氮化物2190条 mean Tc 9.5K, >20K 仅 3.4%; 文献锚点 MoN 33.4K/ZrN 30K
        # 几乎无训练近邻 -> RF 纯外推被压回低Tc均值区(相对低估60%+)。只标注不改数值:
        # 训练集抽样3000条全局 bias=-0.56K(几乎无偏), 全局修正会把正常样本打偏。
        f_up = formula.upper()
        is_metal_nitride = (counts.get('N', 0) > 0) and (counts.get('H', 0) == 0) and (counts.get('O', 0) == 0)
        if is_metal_nitride and tc < 15:
            warnings.append({'id': 'R-20260826-nitride-extrap',
                             'text': '高Tc氮化物外推区: 训练数据集中在该区低Tc(>20K仅3.4%), 预测可能保守偏低, 参考性低'})
        return {'formula': formula, 'counts': counts, 'tc_pred': round(tc, 1),
                'trust_zone': ann['trust_zone'], 'family': ann['family'],
                'warnings': warnings,
                'pressure_boundary': ann['pressure_boundary'],
                'vocab_unknown': ann['vocab_unknown']}, None

def joblib_load(path):
    import joblib
    return joblib.load(path)

ZONE_ICON = {'green': '🟢', 'yellow': '🟡', 'red': '🔴'}

def main():
    args = sys.argv[1:]
    if not args:
        print('用法: python tc_predictor.py "化学式" [--detail]')
        print('示例: python tc_predictor.py "YBa2Cu3O7"')
        print('      python tc_predictor.py "Hg0.7Pb0.3Ba2Ca2Cu3O8" --detail')
        return
    formula = args[0]
    detail = '--detail' in args
    if not os.path.exists(MODEL_PATH):
        print('模型不存在, 请先运行 train_model.py 训练')
        return
    p = TcPredictor()
    res, err = p.predict_formula(formula)
    if err:
        print('错误:', err)
        return
    print('=' * 52)
    print('化学式      : %s' % formula)
    print('组成        : %s' % ', '.join('%s%s' % (k, v) for k, v in sorted(res['counts'].items())))
    print('预测临界温度: %.1f K (摄氏 %.1f°C)' % (res['tc_pred'], res['tc_pred'] - 273.15))
    if res['tc_pred'] >= 77:
        print('判定        : 🟢 高于液氮温度77K, 高温超导潜力')
    elif res['tc_pred'] >= 30:
        print('判定        : 🟡 中低温超导区')
    else:
        print('判定        : 🔴 低温超导区')
    print('信任区      : %s %s (R138 外推校准)' % (ZONE_ICON.get(res['trust_zone'], '?'), res['trust_zone']))
    print('家族        : %s' % res['family'])
    print('压力边界    : %s' % res['pressure_boundary'])
    for w in res['warnings']:
        print('⚠️ [%s] %s' % (w['id'], w['text']))
    print('模型精度    : R2=%.3f RMSE=%.1f K (5折CV)' % (p.meta['cv_r2'], p.meta['cv_rmse']))
    if detail:
        print('-' * 52)
        print('元素明细:')
        for k, v in sorted(res['counts'].items(), key=lambda x: -x[1]):
            print('  %-3s %6.3f' % (k, v))
    print('=' * 52)
    print('⚠️ 预测基于SuperCon 21263种材料数据, 仅供参考, 实际Tc受结构/工艺影响')

if __name__ == '__main__':
    main()
