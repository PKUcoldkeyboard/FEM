# FEM

## 问题描述

根据已知下列非齐次两点边值问题(1.2.28)

$$
\begin{cases}
\boldsymbol{L} u=-\frac{\mathrm{d}}{\mathrm{d} x}\left(p \frac{\mathrm{d} u}{\mathrm{~d} x}\right)+q u=f, a < x < b, \\
u(a)=\alpha, u^{\prime}(b)=\beta,
\end{cases}
$$

与下列变分问题等价：求𝑢 ∈ 𝐻^1, 𝑢(𝑎) = 𝛼，使

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
\text{设}  [a, b]=[-1,1], p(x) \equiv-\left(\pi^{2}-1\right)^{-1}, q(x) \equiv 1, \alpha=0, \beta=-\pi e , \text{以及} \\
$$
$$
f(x)=\frac{2 \pi}{\pi^{2}-1} \cdot \cos (\pi x) \cdot e^{x}
$$

<b>任务1</b>：请认真阅读并完成以下子任务

- 分别取 $ℎ = 0.20, 0.10, 0.05, 0.02$ 时, 将求解域 $[𝑎, 𝑏]$ 等分为长度为$ℎ$的单元或子区间
- 根据上述剖分，就边值问题(1.2.28)和基函数(2.1.16)中 $\varphi_i(𝑥)$ 而设计 $\varphi_0(𝑥)$ ，编程构建相应的 Ritz-Galerkin 方程 (即有限元方程)
- 分别使用<b>高斯消去法</b>和<b>雅可比迭代法</b>(迭代 30 次)，求解上述有限元方程
- 计算得到有限元解 $𝑢_h$ ，绘制有限元解 $𝑢_h$ 的函数图像.

<b>任务2</b>：已知 $u(x)=\sin(\pi x) \cdot e^x$ 是上述边值问题的解析解，针对不同的步长 $h$ 线性方程组解法得到的数值解 $u_h$ ，绘制误差函数 $(u_h − u)$ 的函数图像，且进行观察分析。

## 任务1

### 划分求解域

$[a, b] = [-1, 1], h = 0.2, 0.1, 0.05, 0.02, n = 10, 20, 40, 100$

- $h = 0.2, \quad n = 10,\quad x=[-1, -0.8, \cdots, 0.8, 1]$
  
- $h = 0.1, \quad n = 20, \quad x=[-1, -0.9, \cdots, 0.9, 1]$
  
- $h = 0.05, \quad n = 40, \quad x=[-1, -0.95, \cdots, 0.95, 1]$
  
- $h = 0.02,\quad n = 100, \quad x=[-1, -0.98, \cdots, 0.98, 1]$
  
  ### 分析并构建Ritz-Galerkin方程
  
- 基函数构造，设计 $\varphi_0(x)$ 
  

$$
\varphi_i(x) = \begin{cases}
1+\frac{x-x_i}{h_i}, x_{i-1} \leq x \leq x_i \\
1-\frac{x-x_i}{h_{i+1}}, x_i \leq x \leq x_{i + 1} \\
0, \text{otherwise} \\
\end{cases}
$$

- 左边值条件齐次： $u(a)=\alpha=0$ ，右边值条件非齐次： $u'(b)=\beta$ 
- 则增加基函数：

$$\varphi_0(x) = \begin{cases}
1-\frac{x-x_0}{h_1}, x_0 \leq x \leq x_1\\
0, \text{otherwise} \\
\end{cases}$$

- 满足 $u\in H^1, u(a)=\alpha$ 的试探函数可以写成：

$$u_h(x)= \sum_{i=0}^n \sigma_i \varphi_i(x) = \alpha \varphi_0 (x) + \sum_{i=1}^{n} \sigma_i \varphi_i(x)=\sum_{i=1}^{n} \sigma_i \varphi_i(x), \quad(\sigma_0 = \alpha)$$

- 于是：

$$\begin{aligned}
J(u_h)=&\frac{1}{2} a(u_h, u_h) - (f, u_h) -p(b)\beta(b) u_h(b)\\ 
=&\sum_{i=1}^n \sum_{j=1}^n \sigma_i \sigma_j \frac{a(\varphi_i(x), \varphi_j(x))}{2} -\sum_{i=1}^n \sigma_i (f, \varphi_i(x)) - \sum_{i=1}^n \sigma_i p(b) \beta \varphi_i(b) \\
\end{aligned}$$

- Ritz-Galerkin方程：

$$\begin{aligned}
&\frac{\partial J(u_h)}{\partial \sigma_i}=\sum_{i=1}^n \sigma_i a(\varphi_i, \varphi_j)-(f, \varphi_j)-p(b) \beta \varphi_j(b)=0, \quad j=1, 2, \cdots, n \\
&\Rightarrow \sum_{i=1}^n \sigma_i a(\varphi_i, \varphi_j)=(f, \varphi_j)+p(b) \beta \varphi_j(b), \quad j=1, 2, \cdots, n \\
\end{aligned}$$

- 编程构建Ritz-Galerkin方程并求解，核心代码: `femsolver.py`

### 有限元解的函数图像
- 高斯消去法求解结果
  - Uh-gauss-0.2.png
  - Uh-gauss-0.1.png
  - Uh-gauss-0.05.png
  - Uh-gauss-0.02.png
  - Uh-gauss-0.01.png

- 雅可比迭代法求解结果
  - Uh-jacobi-0.2.png
  - Uh-jacobi-0.1.png
  - Uh-jacobi-0.05.png
  - Uh-jacobi-0.02.png
  - Uh-jacobi-0.01.png

## 任务2：绘制误差函数图像
- 高斯消去法求解结果
  - Error-gauss-0.2.png
  - Error-gauss-0.1.png
  - Error-gauss-0.05.png
  - Error-gauss-0.02.png
  - Error-gauss-0.01.png

- 雅可比迭代法求解结果
  - Error-jacobi-0.2.png
  - Error-jacobi-0.1.png
  - Error-jacobi-0.05.png
  - Error-jacobi-0.02.png
  - Error-jacobi-0.01.png

## 结果分析

- 从高斯消去法求解的结果来看，$u_h$函数近似估计精确解的效果很好，节点处的数值解与精确解的值几乎是重合的，而且随着h的减小误差也不断减少，当h=0.01时，误差的尺度为1e-5至1e-4，基本可以忽略不计。
- 从雅可比迭代法求解的结果来看，$u_h$函数近似估计精确解的效果不太好，节点处的数值解与精确解的值之间误差较大，而且随着h的减少，误差下降到一定程度（1e-2至1e-1）后不再下降。经过程序检验发现造成雅可比迭代不收敛的原因在于对有限元方程构建的总刚度矩阵是一个**非对角占优矩阵**，即不满足雅可比迭代的收敛要求，所以通过雅可比迭代法求解线性方程组 $KU=F$ 无法得到收敛的数值解。

## 许可证
This project is licensed under the MIT License - see the [许可证](许可证) file for details
