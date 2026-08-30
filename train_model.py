# -*- coding: utf-8 -*-
"""超导Tc预测工具 - 模型训练
输入: unique_m.csv (21263条 x 86元素比例 + Tc)
输出: tc_model.joblib (RF模型+元素列顺序+数据统计)
验证: 5折CV R2/RMSE, 保存到工具目录
"""
import csv, os, sys, io, json, time
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score, mean_squared_error
import joblib

root = r'D:\超级黑科技数字人物种'
um = os.path.join(root, 'external_data', 'superconduct_dataset', 'unique_m.csv')
tool_dir = os.path.join(root, 'tc_predictor')
os.makedirs(tool_dir, exist_ok=True)

t0 = time.time()
with open(um, encoding='utf-8') as f:
    rd = csv.reader(f)
    header = next(rd)
    elems = header[:86]
    rows = []
    for r in rd:
        try:
            vals = [float(x) for x in r[:86]]
            tc = float(r[86])
            mat = r[87]
            rows.append((vals, tc, mat))
        except Exception:
            continue
X0 = np.array([x for x, _, _ in rows])
y = np.array([tc for _, tc, _ in rows])
mats = [m for _, _, m in rows]
X = X0 / X0.sum(axis=1, keepdims=True)  # 归一化为比例
print('数据: %d条, Tc min=%.1f max=%.1f mean=%.1f (%.1fs)' % (len(y), y.min(), y.max(), y.mean(), time.time()-t0))

# 5折CV验证
kf = KFold(n_splits=5, shuffle=True, random_state=42)
r2s, rmses = [], []
for tr, te in kf.split(X):
    m = RandomForestRegressor(n_estimators=200, max_depth=20, max_features=0.3,
                              min_samples_leaf=1, n_jobs=-1, random_state=42)
    m.fit(X[tr], y[tr])
    p = m.predict(X[te])
    r2s.append(r2_score(y[te], p))
    rmses.append(np.sqrt(mean_squared_error(y[te], p)))
print('5折CV: R2=%.4f±%.4f RMSE=%.2f±%.2f K' % (np.mean(r2s), np.std(r2s), np.mean(rmses), np.std(rmses)))

# 全量训练最终模型
print('训练全量模型...')
final = RandomForestRegressor(n_estimators=200, max_depth=20, max_features=0.3,
                              min_samples_leaf=1, n_jobs=-1, random_state=42)
final.fit(X, y)
pf = final.predict(X)
print('全量拟合: R2=%.4f RMSE=%.2f K' % (r2_score(y, pf), np.sqrt(mean_squared_error(y, pf))))

# 保存模型+元数据
meta = {
    'elements': elems,
    'n_samples': int(len(y)),
    'tc_min': float(y.min()), 'tc_max': float(y.max()), 'tc_mean': float(y.mean()),
    'cv_r2': float(np.mean(r2s)), 'cv_r2_std': float(np.std(r2s)),
    'cv_rmse': float(np.mean(rmses)),
    'trained_at': time.strftime('%Y-%m-%d %H:%M:%S'),
    'model_type': 'RandomForest(200/20/0.3/1)',
    'note': '输入86维元素比例向量, 输出Tc(K); 比例=化学式原子数归一化',
}
joblib.dump({'model': final, 'meta': meta}, os.path.join(tool_dir, 'tc_model.joblib'))
json.dump(meta, open(os.path.join(tool_dir, 'tc_model_meta.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
print('已保存: %s (%.1f MB)' % (os.path.join(tool_dir, 'tc_model.joblib'), os.path.getsize(os.path.join(tool_dir, 'tc_model.joblib'))/1e6))
