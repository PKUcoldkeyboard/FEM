# FEM
FEM 是一个基于 Taichi 实现的有限元方程求解程序，旨在通过有限元方法解决非齐次两点边值问题。它用于近似求解在给定域内的微分方程解。该项目利用 Taichi 的强大计算能力和可编程特性，将域离散化为较小的单元，并使用适当的基函数构建 Ritz-Galerkin 方程。最后分别通过应用高斯消去法和雅可比迭代等数值技术，求解生成的线性方程组。

## 环境依赖
- Python: 3.9.x
- taichi: 1.6.0

## 源代码
- 单文件: PDE-Project.ipynb

## 结果可视化

## 许可证
This project is licensed under the MIT License. See the LICENSE file for more information.
