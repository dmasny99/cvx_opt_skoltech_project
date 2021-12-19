import numpy as np
import multiprocessing as mp

def z_count(wk, k, df, L, delta, noisy_img):
    wk += ((k + 1) / 2) * df
    zk = -wk / max(L, np.linalg.norm(wk) / delta) + noisy_img
    return wk, zk

def y_count(L, x, noisy_img, df, delta):

    return (L * (x - noisy_img) - df)/ (max(L, np.linalg.norm(L * (x - noisy_img) - df) / delta)) + noisy_img

def gradient_count(df: mp.Value, pobj: mp.Value, i, j, m, x, mu):
    i1 = (i + 1) + j * m
    i2 = i + (j + 1) * m
    i3 = i + j * m
    uij = []
    uij[0] = x[i1] - x[i3]
    uij[1] = x[i2] - x[i3]
    
    c1 = np.linalg.norm(uij)
    pobj += c1
    
    c2 = max(mu, c1)
    uij= uij / c2
        
    df[i1] += uij[0]
    df[i3] -= uij[0]
    df[i2] += uij[1]
    df[i3] -= uij[1]

def test_process(ar, i, j, pobj):
    ar[i, j] = i ** 2 + j ** 2
    pobj += i + j
