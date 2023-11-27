import numpy as np
import matplotlib.pyplot as plt
from functools import partial
from numba import njit, prange, objmode, jit


# setup of 2D q-state Potts model
@njit
def potts_energy(spins, J, h):
    """计算Potts模型的能量

    Args:
        spins: 自旋状态,shape=(N, N)的数组,其中N是格子数目
        J: 相邻自旋互相作用能量强度
        h: 外磁场强度

    Returns:
        系统总能量,标量
    """

    # 获取形状
    N = spins.shape[0]

    # 计算相邻自旋互作用能
    E_J = J * np.sum(spins[:-1, :] == spins[1:, :])
    E_J += J * np.sum(spins[:, :-1] == spins[:, 1:])
    # 计算外磁场能
    E_J += h * np.sum(spins)

    # 求和得到能量
    E = -E_J
    return E


# magnetization
@njit
def potts_magnetization(spins):
    """计算Potts模型的磁化强度

    Args:
        spins: 自旋状态,shape=(N, N)的数组,其中N是格子数目

    Returns:
        系统磁化强度,标量
    """

    # 获取形状
    N = spins.shape[0]

    # 计算磁化强度
    M = np.sum(spins) / N**2
    return M


# Metropolis Monte Carlo Step
@njit(fastmath=True)
def potts_mc_step(spins, J, h, beta, q):
    """一步 Metropolis 蒙特卡洛模拟

    Args:
        spins: 自旋状态,shape=(N, N)的数组,其中N是格子数目
        J: 相邻自旋互相作用能量强度
        h: 外磁场强度
        beta: 逆温度
        q: 自旋状态数目

    Returns:
        spins: 更新后的自旋状态
        rng: 更新后的随机数生成器
    """

    # 获取形状
    N = spins.shape[0]

    # 选取随机格点
    i = np.random.randint(0, N)
    j = np.random.randint(0, N)

    # 建议新的自旋状态
    spins_new = spins.copy()
    spins_new[i, j] = np.random.randint(1, q + 1)

    # 计算能量差
    # dE = potts_energy(spins_new, J, h) - potts_energy(spins, J, h)
    # 直接计算能量差, 由于只改变了一个格点, 只需要计算该格点的能量差
    Eij_old = 0
    Eij_new = 0
    # 计算与上方格点的能量差
    if i > 0:
        if spins[i, j] == spins[i - 1, j]:
            Eij_old -= J
        if spins_new[i, j] == spins_new[i - 1, j]:
            Eij_new -= J
    # 计算与下方格点的能量差
    if i < N - 1:
        if spins[i, j] == spins[i + 1, j]:
            Eij_old -= J
        if spins_new[i, j] == spins_new[i + 1, j]:
            Eij_new -= J
    # 计算与左方格点的能量差
    if j > 0:
        if spins[i, j] == spins[i, j - 1]:
            Eij_old -= J
        if spins_new[i, j] == spins_new[i, j - 1]:
            Eij_new -= J
    # 计算与右方格点的能量差
    if j < N - 1:
        if spins[i, j] == spins[i, j + 1]:
            Eij_old -= J
        if spins_new[i, j] == spins_new[i, j + 1]:
            Eij_new -= J
    # 计算与外磁场的能量差
    Eij_old -= h * spins[i, j]
    Eij_new -= h * spins_new[i, j]
    # 计算能量差
    dE = Eij_new - Eij_old

    # 判断是否接受
    if dE <= 0 or np.random.rand() < np.exp(-beta * dE):
        spins = spins_new

    return spins


# Metropolis Monte Carlo


@jit(forceobj=True)
def potts_mc_mt(T, spins, J, h, q, n_steps, n_step_save=100):
    """Metropolis 蒙特卡洛模拟

    Args:
        T: 温度
        spins: 自旋状态,shape=(N, N)的数组,其中N是格子数目
        J: 相邻自旋互相作用能量强度
        h: 外磁场强度
        n_steps: 模拟步数
        n_step_save: 保存间隔

    Returns:
        spins: 更新后的自旋状态
        rng: 更新后的随机数生成器
    """

    # 进行n_steps步模拟
    spin_history = []
    for i in range(n_steps):
        beta = 1 / (k * T)
        spins = potts_mc_step(spins, J, h, beta, q)
        if i % n_step_save == 0:
            spin_history.append(spins)
    np.save("./Potts_Data/spin_history_{}_{:.2f}.npy".format(q, T), spin_history)
    return spin_history


def potts_simulate_parallel(init_spin, J, h, q, k, T_series, n_steps, n_step_save=100):
    """并行模拟Potts模型并保存轨迹
    Args:
        init_spin: 初始自旋状态,shape=(N, N)的数组,其中N是格子数目
        J: 相邻自旋互相作用能量强度
        h: 外磁场强度
        q: 自旋状态数目
        k: 玻尔兹曼常数
        T_series: 温度序列
        n_steps: 模拟步数
        n_step_save: 保存间隔
    """
    import multiprocessing as mt

    pool = mt.Pool(8)

    # partial function
    potts_mc_partial = partial(
        potts_mc_mt,
        spins=init_spin,
        J=J,
        h=h,
        q=q,
        n_steps=n_steps,
        n_step_save=n_step_save,
    )
    # run simulation
    pool.map(potts_mc_partial, T_series)
    return 0


# perform MC simulation
# setup
q = 3
N = 32
J = 1
h = 0
k = 1
T_series = np.linspace(0.5, 2.5, 100)
n_steps = 5000000
# initial state
spins = np.random.randint(1, q + 1, size=(N, N))
spins = np.array(spins, dtype=np.int8)

# time
import os

if not os.path.exists("./Potts_Data"):
    os.mkdir("./Potts_Data")
potts_simulate_parallel(spins, J, h, q, k, T_series, n_steps, n_step_save=100)
