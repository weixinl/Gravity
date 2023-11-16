#Gravity on a mesh
#plot helper

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.colors as colors
from matplotlib.widgets import Slider


#function to plot 2d z axis slices of the density 
def plot_density(density_field,grid_size=32):
    
    # Set up the figure, the axis, and the plot element
    fig, ax = plt.subplots()

    im = ax.imshow(density_field[:, :, 0], cmap='viridis',norm=colors.LogNorm(vmin=0.001,vmax=np.max(density_field)))
    ax.invert_yaxis()

    #colorbar/density in log scale to better visualize the spread of data
    fig.colorbar(im,ax=ax,label='log(Density)')
    ax.set_title("Density Field Slice (Z = 0)")
    #setting axes ranges
    ax.set_xlim(0, grid_size)
    ax.set_ylim(0, grid_size)
    #setting axes labels
    ax.set_xlabel("X axis")
    ax.set_ylabel("Y axis")

    
    def update(z):
        #can either use the slider to manually go through slices
        idx = int(sliderwave.val) #slider 

        #or automatically have it cycle through all z axis slices
        #idx = int(z)    #auto cycle

        im.set_array(density_field[:, :, idx])
        ax.set_title("Density Field Slice (z = %s)" % int(z))
        
        # if (idx%4==0):
        #     plt.savefig("density_field_slice{}.png".format(np.int32(z)), dpi=1200, bbox_inches="tight")
        return [im] 

    # Sliders
    plt.subplots_adjust(bottom=0.2)
    axwave = plt.axes([0.25, 0.05, 0.5, 0.03])
    
    sliderwave = Slider(axwave, 'Z axis', 0, 31, valinit=0, valfmt='%d')
    sliderwave.on_changed(update)

    #autocycle slices
    #ani = FuncAnimation(fig, update, frames=grid_size)

    plt.show()
    plt.close()

#function to plot particles with given positions
#used for plotting initial distribution or a positions during a timestep
def plot_particles(positions,grid_size=32):
    # plotting particles
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(positions[..., 0], positions[..., 1], positions[..., 2],alpha=0.5, marker='o')
    #setting axes ranges
    ax.set_xlim3d(0, grid_size)
    ax.set_ylim3d(0, grid_size)
    ax.set_zlim3d(0, grid_size)
    #setting axes labels
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    plt.title("Distributing particles with a Gaussian Distribution")

    plt.show()
    #plt.savefig("initial_distribution.png", dpi=1200, bbox_inches="tight")
    plt.close()



