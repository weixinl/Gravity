import numpy as np

def mydft(X):
    '''
    3d fourier transform
    '''
    return

def myinvdft():
    '''
    3d inverse fourier transform
    '''
    return

def fourier(xmat):
    return np.fft.rfftn(xmat)

def inv_fourier(fmat):
    '''
    input: frequency 
    output: time
    '''
    return np.fft.irfftn(fmat)

def solve_poisson(density):
    '''
    solve poisson equation 
    return field
    assuming side of a grid is 1
    '''
    density_freq = fourier(density) # density in frequency 
    N = density.shape[0]
    field_freq = np.zeros((N, N, N))
    for i in range(N):
        for j in range(N):
            for k in range(N):
                field_freq[i, j, k] = - density_freq[i, j, k] / (np.pi * (pow(i + 1, 2) + pow(j + 1, 2) + pow(k + 1, 2)))
    return inv_fourier(field_freq)

if __name__ == "__main__":
    # test
    a = np.zeros((2,2,2))
    a[0,0,0] = 1
    print(solve_poisson(a))

