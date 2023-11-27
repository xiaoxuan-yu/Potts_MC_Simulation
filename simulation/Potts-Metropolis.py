import numpy as np
import matplotlib.pyplot as plt
from functools import partial
from numba import njit, prange, objmode, jit
from tqdm import trange


# Finding neighbors in PBC
@jit
def finding_neighbor(N, i, j):
    return [((i + 1) % N, j), ((1 - i) % N, j), (i, (j + 1) % N), (i, (j - 1) % N)]


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
    """

    # 获取形状
    N = spins.shape[0]

    # 选取随机格点
    i = np.random.randint(0, N)
    j = np.random.randint(0, N)

    # 建议新的自旋状态
    spins_new = spins.copy()
    spins_new[i, j] = np.random.randint(1, q + 1)

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
def potts_mc_mt(T, spins, J, h, q, n_steps, n_step_save=100, path="./Potts_Data"):
    """Metropolis 蒙特卡洛模拟

    Args:
        T: 温度
        spins: 自旋状态,shape=(N, N)的数组,其中N是格子数目
        J: 相邻自旋互相作用能量强度
        h: 外磁场强度
        q: 自旋状态数目
        n_steps: 模拟步数
        n_step_save: 保存间隔
        path: 保存路径

    Returns:
        spin_history: 自旋状态轨迹
    """

    # 进行n_steps步模拟
    spin_history = []
    for i in range(n_steps):
        beta = 1 / (k * T)
        spins = potts_mc_step(spins, J, h, beta, q)
        if i % n_step_save == 0:
            spin_history.append(spins)
    np.save(path + "/spin_history_{}_{:.2f}.npy".format(q, T), spin_history)
    return spin_history


def potts_mc_mh(h, spins, J, T, q, n_steps, n_step_save=100):
    """Metropolis 蒙特卡洛模拟

    Args:
        h: 外磁场强度
        spins: 自旋状态,shape=(N, N)的数组,其中N是格子数目
        J: 相邻自旋互相作用能量强度
        T: 温度
        q: 自旋状态数目
        n_steps: 模拟步数
        n_step_save: 保存间隔

    Returns:
        spin_history: 自旋状态轨迹
    """

    # 进行n_steps步模拟
    spin_history = []
    for i in range(n_steps):
        beta = 1 / (k * T)
        spins = potts_mc_step(spins, J, h, beta, q)
        if i % n_step_save == 0:
            spin_history.append(spins)
    np.save(
        "./Potts_Data_h/spin_history_{}_{}_{:.2f}.npy".format(q, h, T), spin_history
    )
    return spin_history


def potts_simulate_parallel_temp(
    init_spin, J, h, q, k, T_series, n_steps, n_step_save=1000, path="./Potts_Data"
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
        path: 保存路径
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
        path=path,
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

    pool = mt.Pool(6)

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
h_series = np.logspace(-5, 1, 50)
h_series = np.insert(h_series, 0, 0)
n_steps = 10000000
# initial state
spins = np.random.randint(1, q + 1, size=(N, N))
spins = np.array(spins, dtype=np.int8)

# time
import os

# create folder if not exist
if not os.path.exists("./Potts_Data"):
    os.mkdir("./Potts_Data")
if not os.path.exists("./Potts_Data_h"):
    os.mkdir("./Potts_Data_h")
if not os.path.exists("./Potts_Data_critical"):
    os.mkdir("./Potts_Data_critical")

# run simulation
print(" Simulation of different temperatures")
potts_simulate_parallel_temp(spins, J, h, q, k, T_series, n_steps, n_step_save=1000)
print(" Simulation of different magnetic fields")
T = [0.75, 1.2, 2, 3]
for i in trange(len(T)):
    potts_simulate_parallel_h(spins, J, h_series, q, k, T[i], n_steps, n_step_save=1000)
print(" End of simulation.")
