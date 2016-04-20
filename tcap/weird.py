import multiprocessing
import ctypes
import numpy as np

shared_array_base = multiprocessing.Array(ctypes.c_double, 10*10)
shared_array = np.ctypeslib.as_array(shared_array_base.get_obj())
shared_array = shared_array.reshape(10, 10)

# Parallel processing

def my_func(i, shared_array):
    shared_array[i,:] = i

def pool_init(_shared_array, _constans):
	global shared_array, constans
	shared_array = _shared_array
	constans = _constans

def pool_my_func(i):
	my_func(i, shared_array)

if __name__ == '__main__':
    pool = multiprocessing.Pool(8, pool_init, (shared_array, 4))
    for i in np.arange(1000):
        pool.map(pool_my_func, range(10))
    print(shared_array)