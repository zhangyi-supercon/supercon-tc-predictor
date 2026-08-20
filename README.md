# Supercon Tc Predictor 超导临界温度预测工具

AI-powered prediction tool for superconducting critical temperature (Tc) from chemical formula. 基于 21263 条实验数据的超导临界温度智能预测工具。

## 🌐 在线使用（免费试用）

**立即使用：https://tcpredict.top**

- 输入化学式（如 YBa2Cu3O7），数秒得到 Tc 预测值
- 注册即送 10 次免费额度，无需付费即可完整体验
- 充值套餐灵活（30元/10次起），未用次数随时可退

## 📊 模型概况

| 项目 | 数值 |
|------|------|
| 训练数据 | 21,263 条真实实验数据（含 SuperCon 数据库） |
| 交叉验证 R² | 0.922 |
| 支持体系 | 铜基、铁基、氢化物、传统超导体等 |
| 输出 | 临界温度 Tc（K）+ 置信参考 |

## 🔬 适用场景

- 材料科研人员：快速筛选候选材料，缩小实验范围
- 课题组：做初步预判，节省第一性原理计算成本
- 超导爱好者：探索新材料体系的 Tc 规律

## ⚙️ API 用法

注册后获取 API Key（X-API-Key 请求头）：

```bash
curl "https://tcpredict.top/api/predict?formula=YBa2Cu3O7" \
  -H "X-API-Key: YOUR_API_KEY"
```

```python
import requests
url = "https://tcpredict.top/api/predict"
r = requests.get(url, params={"formula": "YBa2Cu3O7"},
                  headers={"X-API-Key": "YOUR_API_KEY"})
print(r.json()["tc"])
```

## ⚠️ 说明

- 预测结果仅供科研参考，不能替代第一性原理计算与实验室验证
- 数据来源：公开超导数据库（SuperCon 等）与文献实验数据
- 持续更新中

## 📮 联系

- 官网：https://tcpredict.top
- 注册邮箱：382776397@qq.com

---

*Made with ❤️ for the superconductivity research community*
