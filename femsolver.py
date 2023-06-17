import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.integrate import quad


class FEMSolver():
    def __init__(self, h=None, solver='gauss'):
        self.h = h

        self.a, self.b = -1.0, 1.0
        self.alpha, self.beta = 0.0, -np.pi * np.exp(1)
        self.n = int((self.b - self.a) / self.h)

        self.X = np.linspace(self.a, self.b, self.n + 1)
        self.K = np.zeros((self.n + 1, self.n + 1))
        self.F = np.zeros(self.n + 1)
        self.U = np.zeros(self.n + 1)
        self.solver = solver

    def f(self, x):
        return (2 * np.pi / (np.pi**2 - 1)) * np.cos(np.pi * x) * np.exp(x)

    def p(self, x):
        return -(np.pi**2 - 1)**-1

    def q(self, x):
        return 1.0

    def exact(self, x):
        return np.sin(np.pi * x) * np.exp(x)

    def phi(self, i, x):
        if x >= self.X[i - 1] and x <= self.X[i]:
            return 1 + (x - self.X[i]) / self.h
        elif x >= self.X[i] and x <= self.X[i + 1]:
            return 1 - (x - self.X[i]) / self.h
        else:
            return 0

    def dphi(self, i, x):
        if x >= self.X[i - 1] and x <= self.X[i]:
            return 1 / self.h
        elif x >= self.X[i] and x <= self.X[i + 1]:
            return -1 / self.h
        else:
            return 0

    def aij(self, i, j):
        term1 = quad(lambda x: self.p(x) * self.dphi(i, x) * self.dphi(j, x),
                     max(self.a, self.X[i]-self.h), min(self.b, self.X[i]+self.h))[0]
        term2 = quad(lambda x: self.q(x) * self.phi(i, x) * self.phi(j, x),
                     max(self.a, self.X[i]-self.h), min(self.b, self.X[i]+self.h))[0]
        return term1 + term2

    # 总刚度矩阵
    def build_K(self):
        for i in range(self.n + 1):
            for j in range(self.n + 1):
                self.K[i, j] = self.aij(i, j)
        # 左边值为0
        self.K = self.K[1:, 1:]

    # 总载荷向量
    def build_F(self):
        for i in range(self.n + 1):
            self.F[i] = quad(lambda x: self.f(x)*self.phi(i, x), max(self.a, self.X[i]-self.h),
                             min(self.b, self.X[i]+self.h))[0] + self.p(self.b)*self.phi(i, self.b)*self.beta
        # 左边值为0
        self.F = self.F[1:]

    def solve(self):
        self.build_K()
        self.build_F()
        self.U = np.zeros(self.n)

        if self.solver == 'gauss':
            self.U = self.gauss_solve(self.K, self.F)
        elif self.solver == 'jacobi':
            self.U = self.jacobi_solve(self.K, self.F)

        # 增加边界条件u(a)=0
        self.U = np.concatenate(([0], self.U))

    # 高斯消去法
    def gauss_solve(self, K, F):
        n = len(K)
        nf = F.shape[0]
        # 1. 消元
        for i in range(0, n - 1):
            for j in range(i + 1, n):
                lam = K[j, i] / K[i, i]
                K[j, i:n] = K[j, i:n] - lam * K[i, i:n]
                F[j] = F[j] - lam * F[i]

        # 2. 回代
        for i in range(n - 1, -1, -1):
            F[i] = (F[i] - np.dot(K[i, i + 1:n], F[i + 1:n])) / K[i, i]

        return F

    # 雅可比迭代法
    def jacobi_solve(self, K, F, max_iter=30):
        # 首先判断是强对角占优矩阵、弱对角占优矩阵还是不对角占优矩阵
        if np.all(np.abs(np.diag(K)) > np.sum(np.abs(K), axis=1) - np.abs(np.diag(K))):
            print('>>> 强对角占优矩阵')
        elif np.all(np.abs(np.diag(K)) >= np.sum(np.abs(K), axis=1) - np.abs(np.diag(K))):
            print('>>> 弱对角占优矩阵')
        else:
            print('>>> 非对角占优矩阵')
            
        L = np.array(np.tril(K, -1))
        U = np.array(np.triu(K, 1))
        D_inv = np.diag(1 / np.diag(K))
        x = np.zeros(self.n)

        for _ in range(max_iter):
            x_new = x
            x = D_inv.dot(F - L.dot(x) - U.dot(x))
            if abs(x - x_new).max() < 1e-6:
                break
        return x.flatten()

    def plot(self):
        sns.set_context('paper')
        plt.figure(figsize=(10, 6), dpi=300)
        # 画出U的图像
        sns.lineplot(x=self.X, y=self.U, label='$u_h$', color='green')
        plt.legend()

        plt.xlabel('x')
        plt.ylabel('$u_h(x)$')
        plt.title('FEM Solution')
        plt.savefig(f'Uh-{self.solver}-{self.h}.png')

        # 清空plt，画误差图像
        errors = np.abs(self.U - self.exact(self.X))
        plt.clf()
        plt.figure(figsize=(10, 6), dpi=300)
        sns.lineplot(x=self.X, y=errors, label='Error', color='green')
        plt.legend()
        plt.xlabel('x')
        plt.ylabel('Error')
        plt.title('Error between FEM Solution and Exact Solution')
        plt.savefig(f'Error-{self.solver}-{self.h}.png')
