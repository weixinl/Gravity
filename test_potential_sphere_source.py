#!/usr/bin/env python
# coding: utf-8

# In[10]:


import numpy as np
import densities 
import plotter
import matplotlib.pyplot as plt


# In[16]:


def test_potential_sphere(center=(16,16,16),n_particles=500,gridsize=32):
    positions = densities.distribute_spherical(center,radius=5,grid_size=gridsize, num_particles=n_particles)
    plotter.plot_particles(positions,save_fig=True,title="test potential from spherical source",
                           filename="imgs/test_potential_sphere.png",particle_size=10)
    
    density = cic_density(positions, grid_size=32,cell_size =1,mass=1.)
    plotter.plot_density_slice(density,time=2,axis='z',slice_position=16,grid_size=32,save_plot=False)
     
        
        
def test_potential_abs_error(center=(16,16,16),gridsize=32):
    X,Y,Z = np.mgrid[0:boxsize,0:boxsize,0:boxsize]
    
    posx = center[0]
    posy = center[1]
    posz = center[2]
    eps = 1e-17
    r = np.sqrt((X-posx)**2+(Y-posy)**2+(Z-posz)**2)   

    rho = np.sum(density)

    ana_V = -density/(r)+eps
    ana_V[posx,posy,posz]=np.NaN

    error = abs(ana_V[:,:,16]-potential[:32,:32,16])
    
    fig, ax = plt.subplots(figsize=(8,6))
    pic2 = plt.imshow(error)
    fig.colorbar(pic2)
    plt.xlabel('Distance from x=0 (m)',fontsize=16)
    plt.ylabel('Distance from y=0 (m)',fontsize=16)
    plt.title('Absolute error',fontsize=20)

if __name__=="__main__":
    test_potential_sphere()
    test_potential_abs_error()


# In[15]:




