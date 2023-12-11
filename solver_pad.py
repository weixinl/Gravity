#!/usr/bin/env python
# coding: utf-8

# In[ ]:


## old

#Gravity on a mesh
import matplotlib.pyplot as plt
import numpy as np

from densities import *
from distribute import *
from plotter import *
from integrator import *

def main():
    # Generate particles and density field
    grid_size = 32
    num_particles = 3
    center = (grid_size/2, grid_size/2,grid_size/2)
    a = 5
    b_to_a = 0.5
    c_to_a = 0.8

    #particles = distribute_spherical(radius=10,num_particles=num_particles)
    particles = distribute_particles(center, a, b_to_a, c_to_a,grid_size=grid_size,num_particles=num_particles)
    # particles = np.zeros(np.shape(particles))
    # particles[0] = [14,16,16]
    # particles[1] = [18,16,16]
    # particles[2] = [16,16,20]
    # particles[3] = [16,15,15]

    density = cic_density(particles,grid_size=grid_size)
    #plot_particles(particles)
    #plot_density(density,grid_size)

    #array to store all positions at each timestep
    moves = np.array(particles)
    #adding a new axis to store timesteps: moves[0] gives positions of all particles at t=0
    moves = moves[np.newaxis,...]

    #testing out with zeroes first
    velocity = np.zeros(np.shape(particles))

    tsteps = 100

    for i in range(tsteps):
        particles,velocity= integrate(particles,velocity,density,time_step=0.5)
        density = cic_density(particles,grid_size=grid_size)
        # enforcing toroidal boundary conditions
        particles = particles%32
        moves = append_new_array(moves,particles)

    plot_particle_motion(moves, tsteps, interval = 3.2,save_gif=False,name="3particles.gif",particle_size=10)
    #plot_motion_from_save("test.txt",num_particles)

if __name__ == '__main__':
    main()


# In[75]:


# Peak is at corners
def green_func_pad2D(N):
    x = np.linspace(0,N,N,endpoint=False)
    y = np.linspace(0,N,N,endpoint=False)
    
    X,Y = np.meshgrid(x,y,indexing='ij')
    eps = 1e-17
    r = np.sqrt(X**2+Y**2)
    eps = 1e-10
    g = np.where(r!=0,1/(r+eps),0)
    g[0,0]=1
    
    
    for i in range(N):
        for j in range(N):
            g[i,j] = g[min(i,N-1-i),min(j,N-1-j)]

    return g


# In[ ]:


def solve_poisson_g(density,gx):
    '''
    solve poisson equation 
    return potential
    assuming side of a grid is 1
    '''
    density_freq = fourier(density) # density in frequency 

    N = density.shape[0]
#     plt.imshow(np.real(density_freq))
    potential_freq = density_freq * gx
#     potential_freq = density_freq * gx
    return inv_fourier(potential_freq)

