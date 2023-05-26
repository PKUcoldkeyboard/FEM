# FEM

FEM 是一个基于 Taichi 实现的有限元方程求解程序，旨在通过有限元方法解决非齐次两点边值问题。它用于近似求解在给定域内的微分方程解。该项目利用 Taichi 的强大计算能力和可编程特性，将域离散化为较小的单元，并使用适当的基函数构建 Ritz-Galerkin 方程。最后分别通过应用高斯消去法和雅可比迭代等数值技术，求解生成的线性方程组。

## 环境依赖

- Python: 3.9.x
- taichi: 1.6.0

## 源代码

- 单文件: PDE-Project.ipynb、

## 问题描述

根据已知下列非齐次两点边值问题(1.2.28)

$$
\left\{\begin{array}{l}

\boldsymbol{L} u=-\frac{\mathrm{d}}{\mathrm{d} x}\left(p \frac{\mathrm{d} u}{\mathrm{~d} x}\right)+q u=f, a<x<b, \\

u(a)=\alpha, u^{\prime}(b)=\beta,

\end{array}\right.

$$

与下列变分问题等价：求$u_{*} \in H^1, u(a)=\alpha$，使

$$
J(u_{*})=\min\limits_{u \in H^1 \atop u(a) = \alpha} J(u)

$$

其中

$$
\begin{array}{c}

J(u)=\frac{1}{2} a(u, u)-(f, u)-p(b) \beta u(b)

\end{array}

$$

$$
a(u, v)=\int_{a}^{b}\left(p \frac{d u}{d x} \frac{d v}{d x}+q u v\right) d x,\quad(g, u)=\int_{a}^{b} g u d x . \\

$$

$$
\text{设}  [a, b]=[-1,1], p(x) \equiv-\left(\pi^{2}-1\right)^{-1}, q(x) \equiv 1, \alpha=0, \beta=-\pi e , \text{以及} \\

$$

$$
f(x)=\frac{2 \pi}{\pi^{2}-1} \cdot \cos (\pi x) \cdot e^{x}

$$

<b>任务1</b>：请认真阅读并完成以下子任务

- 分别取$ℎ = 0.20, 0.10, 0.05, 0.02$时, 将求解域$[𝑎, 𝑏]$等分为长度为$ℎ$的单元或子区间

- 根据上述剖分，就边值问题(1.2.28)和基函数(2.1.16)中$\varphi_i(𝑥)$而设计$\varphi_0(𝑥)$，编程构建相应的 Ritz-Galerkin 方程 (即有限元方程)

- 分别使用<b>高斯消去法</b>和<b>雅可比迭代法</b>(迭代 30 次)，求解上述有限元方程

- 计算得到有限元解$𝑢_h$，绘制有限元解$𝑢_h$的函数图像.

<b>任务2</b>：已知$u(x)=\sin(\pi x) \cdot e^x$是上述边值问题的解析解，针对不同的步长$h$线性方程组解法得到的数值解$u_h$，绘制误差函数$(u_h − u)$的函数图像，且进行观察分析。

## 结果可视化



## 许可证

This project is licensed under the MIT License. See the LICENSE file for more information.
