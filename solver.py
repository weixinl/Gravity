import numpy as np
epsilon = 1e-7 #for numerical stability

def fourier(xmat):
    return np.fft.fftn(xmat)

def inv_fourier(fmat):
    '''
    input: frequency 
    output: time
    '''
    return np.real(np.fft.ifftn(fmat))

def solve_poisson_old(density):
    '''
    solve poisson equation 
    return potential
    assuming side of a grid is 1
    '''
    density_freq = fourier(density) # density in frequency 
    N = density.shape[0]
    potential_freq = np.zeros((N, N, N),dtype=np.cfloat)
    for i in range(N):
        for j in range(N):
            for k in range(N):
                potential_freq[i, j, k] = density_freq[i, j, k] * 2 * np.pi/ (np.cos(2 * np.pi * i / N) + np.cos(2 * np.pi * j / N) + np.cos( 2 * np.pi * k / N) - 3 + epsilon)
    return inv_fourier(potential_freq)
    
def green_func_pad3D(padded_boxlength):
    '''
    Padded Green's function that eleminiates periodic effects of FT
    :param boxlength: np.array 
        Length of the simulating box
    '''
    # Create a larger box that is doubled size
    N = padded_boxlength
    x = np.linspace(0,N,N,endpoint=False)
    y = np.linspace(0,N,N,endpoint=False)
    z = np.linspace(0,N,N,endpoint=False)
    
    # Create a 3D meshgrid
    X,Y,Z = np.meshgrid(x,y,z,indexing='ij')
    
    # Create 1/r function throughout the space
    epsi = 1e-17
    r = np.sqrt(X**2+Y**2+Z**2)
    g = np.where(r!=0,1/(r+epsi),0)
    # Green's function = 1 at the origin
    g[0,0,0]=1
    
    # create mirrored Green's function for the purpose of canclling out FT periodic effect
    for i in range(N):
        for j in range(N):
            for k in range(N):
                g[i,j,k] = g[min(i,N-i),min(j,N-j),min(k,N-k)]
    return g

def patch_density_box(density):
    '''
    Take in the original density with size NxNxN. Patch the density to 2N*2N*2N box and filled the additional space with 0 mass density. 
    The original box is placed at the origin such that Patched_box[0:N,0:N,0:N] is the original box. 
    
    :param density: np.array
        array of shape (grid_size,grid_size,grid_size) with density field due to all the particles
    
    :return density: np.array
            returns array of shape (2*grid_size,2*grid_size,2*grid_size) with density field due to all the particles
    '''
    N = density.shape[0]
    patched_density = np.zeros((2*N,2*N,2*N))
    patched_density[0:N,0:N,0:N] = density
    return patched_density

def solve_poisson(density):
    '''
    solve poisson equation 
    return potential
    assuming side of a grid is 1
    '''
    '''
    solve poisson equation 
    return potential
    assuming side of a grid is 1
    '''
    # Extract original box length 
    orig_N = density.shape[0]
    
    # Create a patched mass density grid
    patched_density = patch_density_box(density)
    
    # Extract padded box length 
    N = patched_density.shape[0]
    
    # Create a Green's function of the padded box size
    g = green_func_pad3D(N)
    
    # Green's function in Fourier space
    g_freq = fourier(g)
    
    # density in Fourier space
    density_freq = fourier(patched_density) 
    
    # Initialize potential (with padded box size) in Fourier space
    potential_freq = np.zeros((N, N,N),dtype=np.cfloat)

    # Calculate the potential by convoluting density_freq and g_freq
    potential_freq = density_freq* g_freq
                
    # Take inverse FT to find the potential in coordinate space
    potential_pad = -inv_fourier(potential_freq)
    
    # Extract the potential with the original size
    potential = potential_pad[0:orig_N,0:orig_N,0:orig_N]
    return potential
    
if __name__ == "__main__":
    # test
    a = np.zeros((2,2,2))
    a[0,0,0] = 1
    print(solve_poisson(a))

