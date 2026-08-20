#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
超导Tc预测器 示例代码
Supercon Tc Predictor - Example Usage

在线使用: https://tcpredict.top
注册后即可获得免费 10 次试用额度与专属 API Key。
"""

import requests

API_URL = "https://tcpredict.top/api/predict"
API_KEY = "YOUR_API_KEY"  # 在 https://tcpredict.top/account 查看（注册邮箱可找回）


def predict(formula: str, api_key: str = API_KEY) -> dict:
    """输入化学式，返回超导临界温度预测结果。"""
    resp = requests.get(
        API_URL,
        params={"formula": formula},
        headers={"X-API-Key": api_key},
        timeout=30,
    )
    resp.raise_for_status()
    return resp.json()


if __name__ == "__main__":
    # 示例 1：钇钡铜氧（第一个 Tc 超过液氮沸点的高温超导体）
    r = predict("YBa2Cu3O7")
    print("YBa2Cu3O7 ->", r)

    # 示例 2：铜基超导体 La2CuO4
    r = predict("La2CuO4")
    print("La2CuO4  ->", r)

    # 示例 3：铁基超导体
    r = predict("BaFe2As2")
    print("BaFe2As2 ->", r)
