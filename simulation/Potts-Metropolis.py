import numpy as np
import matplotlib.pyplot as plt
from functools import partial
from numba import njit, prange, objmode, jit


# Finding neighbors in PBC
@jit
def finding_neighbor(N, i, j):
    return [((i + 1) % N, j), ((1 - i) % N, j), (i, (j + 1) % N), (i, (j - 1) % N)]


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
    E_J = J * np.sum(spins == np.roll(spins, 1, 0))
    E_J += J * np.sum(spins == np.roll(spins, 1, 1))
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
    # 计算与相邻格点的能量差
    neighbors = finding_neighbor(N, i, j)
    for neighbor in neighbors:
        Eij_old -= J * (spins[i, j] == spins[neighbor[0], neighbor[1]])
        Eij_new -= J * (spins_new[i, j] == spins[neighbor[0], neighbor[1]])

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


def potts_mc_mh(h, spins, J, T, q, n_steps, n_step_save=100):
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
    np.save(
        "./Potts_Data_h/spin_history_{}_{:.2f}_{:.2f}.npy".format(q, h, T), spin_history
    )
    return spin_history


def potts_simulate_parallel_temp(
    init_spin, J, h, q, k, T_series, n_steps, n_step_save=1000
):
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


def potts_simulate_parallel_h(
    init_spin, J, h_series, q, k, T, n_steps, n_step_save=100
):
    """并行模拟Potts模型并保存轨迹
    Args:
        init_spin: 初始自旋状态,shape=(N, N)的数组,其中N是格子数目
        J: 相邻自旋互相作用能量强度
        h_series: 外磁场强度序列
        q: 自旋状态数目
        k: 玻尔兹曼常数
        T: 温度
        n_steps: 模拟步数
        n_step_save: 保存间隔
    """
    import multiprocessing as mt

    pool = mt.Pool(8)

    # partial function
    potts_mc_partial = partial(
        potts_mc_mh,
        spins=init_spin,
        J=J,
        T=T,
        q=q,
        n_steps=n_steps,
        n_step_save=n_step_save,
    )
    # run simulation
    pool.map(potts_mc_partial, h_series)
    return 0


# perform MC simulation
# setup
q = 3
N = 32
J = 1
h = 0
k = 1
T_series = np.linspace(0.5, 2.5, 100)
h_series = np.logspace(0, 1, 20)
n_steps = 10000000
# initial state
spins = np.random.randint(1, q + 1, size=(N, N))
spins = np.array(spins, dtype=np.int8)

# time
import os

if not os.path.exists("./Potts_Data"):
    os.mkdir("./Potts_Data")
if not os.path.exists("./Potts_Data_h"):
    os.mkdir("./Potts_Data_h")
print(" Simulation of different temperatures")
potts_simulate_parallel_temp(spins, J, h, q, k, T_series, n_steps, n_step_save=1000)
print(" Simulation of different magnetic fields")
potts_simulate_parallel_h(spins, J, h_series, q, k, 0.75, n_steps, n_step_save=1000)
potts_simulate_parallel_h(spins, J, h_series, q, k, 1.2, n_steps, n_step_save=1000)
potts_simulate_parallel_h(spins, J, h_series, q, k, 2, n_steps, n_step_save=1000)
potts_simulate_parallel_h(spins, J, h_series, q, k, 3, n_steps, n_step_save=1000)
print(" Simulation finished")
