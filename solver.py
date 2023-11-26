import numpy as np
epsilon = 1e-7 #for numerical stability

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
    return np.fft.fftn(xmat)

def inv_fourier(fmat):
    '''
    input: frequency 
    output: time
    '''
    return np.fft.irfftn(fmat)

def solve_poisson(density):
    '''
    solve poisson equation 
    return potential
    assuming side of a grid is 1
    '''
    density_freq = fourier(density) # density in frequency 
    N = density.shape[0]
    potential_freq = np.zeros((N, N, N))
    for i in range(N):
        for j in range(N):
            for k in range(N):
                potential_freq[i, j, k] = density_freq[i, j, k] * 2 * np.pi/ (np.cos(2 * np.pi * i / N) + np.cos(2 * np.pi * j / N) \
                                                                          + np.cos( 2 * np.pi * k / N) - 3 + epsilon)
    return inv_fourier(potential_freq)

if __name__ == "__main__":
    # test
    a = np.zeros((2,2,2))
    a[0,0,0] = 1
    print(solve_poisson(a))

