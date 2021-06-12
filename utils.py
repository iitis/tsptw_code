import numpy as np
import math

def write_npz_file(dir, filename, **kwargs):
    file = f"{dir}/{filename}.npz"
    np.savez(file, **kwargs)


def load_npz_file(dir, filename, data):
    file = f"{dir}/{filename}.npz"
    return np.load(file, allow_pickle='TRUE')[data]


def smallest(array, n):
    array = np.array(array)
    array = array.flatten()

    return sorted(array)[n - 1]


def find_bound(Q, t, arr):
    d = arr.max()
    for i in range(1, t):
        d -= smallest(Q, len(Q) + i)
    if d <= 0:
        return 1
    return math.floor(math.log(d, 2)) + 1


def find_bound_u(Q, t, arr):
    d = arr.max()
    for i in range(1, t):
        d -= smallest(Q, len(Q) + i)
    if d <= 0:
        return 1
    return int(d + 1)
