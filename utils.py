import numpy as np
import math

def fibonacci_sphere(samples=1000):

    points = []
    phi = math.pi * (math.sqrt(5.) - 1.)  # golden angle in radians

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = math.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = math.cos(theta) * radius
        z = math.sin(theta) * radius

        points.append((x, y, z))

    return points

def calc_array_len(arr):
    '''
    calculate the length of an array
    '''
    return np.sqrt(np.sum(arr * arr))

def normalize_array(arr):
    '''
    :param direction: tuple
    '''
    arr_len = calc_array_len(arr)
    assert arr_len > 1e-8 # length is not 0
    return arr/arr_len

def rotate_x(arr, theta = np.pi/2):
    '''
    rotate about x axis
    :param arr: np.array
            array to rotate 
    :param theta: float
            rotation angle
    '''
    rx = np.zeros((3, 3))
    rx[0, 0] = 1
    rx[1, 1] = np.cos(theta)
    rx[1, 2] = -np.sin(theta)
    rx[2, 1] = np.sin(theta)
    rx[2, 2] = np.cos(theta)
    return np.dot(rx, arr)

def rotate_y(arr, theta = np.pi/2):
    '''
    rotate about y axis
    :param arr: np.array
            array to rotate 
    :param theta: float
            rotation angle
    '''
    ry = np.zeros((3, 3))
    ry[0, 0] = np.cos(theta)
    ry[0, 2] = np.sin(theta)
    ry[1, 1] = 1
    ry[2, 0] = -np.sin(theta)
    ry[2, 2] = np.cos(theta)
    return np.dot(ry, arr)

def rotate_by_axis(arr, axis, theta):
    '''
    rotate about a specified axis
    :param arr: np.array
            array to rotate 
    :param axis: np.array
            rotation axis
    :param theta: float
            rotation angle
    '''
    
    # normalize axis arr
    axis = normalize_array(axis)
    
    r = np.zeros((3, 3))
    ux = axis[0] 
    uy = axis[1]
    uz = axis[2]
    cosval = np.cos(theta)
    sinval = np.sin(theta)
    r[0, 0] = cosval + ux * ux * (1 - cosval)
    r[0, 1] = ux * uy * (1 - cosval) - uz * sinval
    r[0, 2] = ux * uz * (1 - cosval) + uy * sinval
    r[1, 0] = uy * ux * (1 - cosval) + uz * sinval
    r[1, 1] = cosval + uy * uy * (1 - cosval)
    r[1, 2] = uy * uz * (1 - cosval) - ux * sinval
    r[2, 0] = uz * ux * (1 - cosval) - uy * sinval
    r[2, 1] = uz * uy * (1 - cosval) + ux * sinval
    r[2, 2] = cosval  + uz * uz * (1 - cosval)

    return np.dot(r, arr)

if __name__ == "__main__":
    # test
    arr = np.zeros(3)
    arr[0] = 1
    print(f"old array: {arr}") # [1, 0, 0]
    arr_new = rotate_by_axis(arr, np.array([0, 1, 0]), np.pi/2)
    print(f"new array: {arr_new}") # should be [0, 0, -1]