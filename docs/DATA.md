# 数据说明 (Data Description)

## 训练数据

- **总量**: 21,263 条真实实验记录
- **来源**: SuperCon 公开超导数据库 + 已发表文献中的实验数据
- **字段**: 化学式（元素组成）、合成/测量条件、实验临界温度 Tc

## 主要覆盖体系

| 体系 | 示例化学式 | 说明 |
|------|-----------|------|
| 铜基超导体 | YBa2Cu3O7, La2CuO4, Bi2Sr2CaCu2O8 | 高温超导主力体系 |
| 铁基超导体 | BaFe2As2, SmFeAsO | 2008 年后兴起 |
| 氢化物 | LaH10, H3S | 高压超导前沿 |
| 传统超导体 | Nb3Sn, MgB2 | BCS 理论体系 |

## 模型指标

- **交叉验证 R²**: 0.922
- **特征**: 元素组成 + 结构相关信息
- **方法**: 机器学习回归模型（详见官网更新说明）

## 使用注意

1. 预测结果仅供科研参考与候选筛选，**不能替代第一性原理计算与实验室验证**
2. 部分体系训练样本较少，预测置信度较低（API 返回中会给出参考信息）
3. 数据持续更新中

## 建议引用

使用本工具辅助研究时，建议在论文/报告中注明：

> Tc prediction performed with Supercon Tc Predictor (https://tcpredict.top), trained on SuperCon database and literature experimental data.

---

更多信息: https://tcpredict.top
