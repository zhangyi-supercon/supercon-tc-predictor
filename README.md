# Supercon Tc Predictor

**AI-powered prediction of superconducting critical temperature (Tc) from chemical formula.**
基于 21,263 条真实实验数据的超导临界温度预测工具 — 分布内极准，外推诚实标注。

![predict vs experiment](predict_vs_experiment.png)

## Why this project exists

Predicting Tc from a chemical formula is hard: synthesis is expensive, and theory-based screening
doesn't always transfer. This project trains a Random Forest on 21,263 real experimental records
(from SuperCon and literature) to answer one question fast:

> Given a chemical formula, what Tc would experiments likely measure — and **how much should we trust the answer?**

## Model summary

| Metric | Value |
|---|---|
| Training data | 21,263 experimental records (SuperCon + literature) |
| Algorithm | Random Forest (200 trees) |
| Cross-validation R² | **0.922** |
| RMSE | **9.6 K** |
| Feature space | 81 chemical descriptors (composition statistics: mean / weighted mean / geometric mean / entropy / range / std) |
| Interface | chemical formula → Tc (K) |

## Signature cases (real predictions)

| Formula | Predicted | Experimental | Note |
|---|---|---|---|
| HgBa₂Ca₂Cu₃O₈ (Hg-1223) | 126.9 K | ~135 K | within ~6% of the ambient-pressure record holder |
| YBa₂Cu₃O₇ (YBCO) | 85.5 K | ~92 K | within 7 K |
| MgB₂ | 35.7 K | 39 K | within 3.3 K |
| Nb₃Sn | 17.2 K | 18 K | within 1 K |
| MoN | 7.7 K | 33.4 K (claimed) | **out-of-distribution: honest miss** |
| LaH₁₀ | 64.1 K | ~250 K (high pressure) | **flagged: extrapolation only, lower-bound reference** |

The last two rows are the point: the model is highly accurate **inside** its training distribution,
and it **knows** when it is outside it. Instead of overconfident numbers, it emits explicit warnings
on extrapolation regions (see `boundary_rules.json`).

## Quick start (local)

```bash
pip install numpy pandas scikit-learn joblib matplotlib

# predict a single formula
python tc_predictor.py "YBa2Cu3O7"
python tc_predictor.py "HgBa2Ca2Cu3O8" --detail
```

### Example output

```
化学式      : HgBa2Ca2Cu3O8
预测临界温度: 126.9 K (摄氏 -146.2°C)
判定        : 🟢 高于液氮温度77K, 高温超导潜力
信任区      : 🔴 red (外推校准)
家族        : cuprate_Hg
模型精度    : R2=0.922 RMSE=9.6 K (5折CV)
⚠️ [B1] 预测>=100K 的输出必须标记 '外推区, 需文献验证'
⚠️ [B8] 高Tc端统计模型可能系统性低估, 预测应视为下界参考
```

## Repository layout

```
├── tc_predictor.py          # main CLI: formula → Tc + boundary warnings
├── boundary_rules.json      # extrapolation / family / pressure boundary rules
├── train_model.py           # training pipeline (RandomForest, 5-fold CV)
├── example.py               # online API usage example
├── docs/
│   ├── DATA.md              # dataset description & preprocessing
│   └── tc-prediction-practice.md  # practical notes from real usage
└── predict_vs_experiment.png
```

## Online demo

Try it live (free trial, no credit card): **https://tcpredict.top**

- Input a formula (e.g. `YBa2Cu3O7`), get Tc + confidence reference in seconds
- 10 free credits on registration; unused credits refundable

## Honesty policy

- Predictions are research references only — they do **not** replace first-principles calculations or lab verification.
- Predictions ≥100 K are always flagged as "extrapolation region, needs literature verification."
- High-Tc extrapolations are treated as **lower-bound references** (statistical models systematically under-estimate in the high-Tc regime — see Xie et al., npj Comput. Mater. 2022).
- The model never sorts candidates by raw predicted Tc in the ≥100 K region (known anti-correlation, B4 rule).

## Contact

- Website: https://tcpredict.top
- Email: 382776397@qq.com

---

*Made with ❤️ for the superconductivity research community*
