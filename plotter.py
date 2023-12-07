#Gravity on a mesh
#plot helper

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import matplotlib.colors as colors
from matplotlib.widgets import Slider
from densities import *

#function to plot particles with given positions
#used for plotting initial distribution or a positions during a timestep
def plot_particles(positions,grid_size=32,save_fig=False,particle_size=0.1):
    # plotting particles
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter3D(positions[..., 0], positions[..., 1], positions[..., 2], s=particle_size)

    #setting axes ranges
    ax.set_xlim3d(0, grid_size)
    ax.set_ylim3d(0, grid_size)
    ax.set_zlim3d(0, grid_size)

    #setting axes labels
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_zlabel('Z-axis')
    #making axis labels smaller
    ax.xaxis.label.set_size(7)
    ax.yaxis.label.set_size(7)
    ax.zaxis.label.set_size(7)
    #making axis ticks smaller
    ax.tick_params(axis='both', which='major', labelsize=5)
    plt.title("Distributing particles with a Gaussian Distribution")

    #making the axis labels closer to the axis
    ax.xaxis.labelpad = -10
    ax.yaxis.labelpad = -10
    ax.zaxis.labelpad = -10

    #making the axis ticks closer to the axis
    ax.tick_params(axis='x', pad=-5)
    ax.tick_params(axis='y', pad=-5)
    ax.tick_params(axis='z', pad=-4)

    ax.scatter3D(positions[..., 0], positions[..., 1], 0 , c='grey', s=0.01)
    ax.scatter3D(0, positions[..., 1], positions[..., 2], c='grey', s=0.01)
    ax.scatter3D(positions[..., 0], grid_size , positions[..., 2], c='grey', s=0.01)

    if save_fig:
        plt.savefig("initial_distribution.png", dpi=1200, bbox_inches="tight")
    else:
        plt.show()

    plt.close()


#function to plot z axis slices of the density field
def plot_density(density_field,grid_size=32,slider=False,save_slices=False,interval=4):
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
        if slider:
            idx = int(sliderwave.val)
        else:
            idx = int(z)
        im.set_array(density_field[:, :, idx])
        ax.set_title("Density Field Slice (z = %s)" % int(z))
        if save_slices:
            if (idx%interval==0):
                plt.savefig("density_field_slice{}.png".format(np.int32(z)), dpi=1200, bbox_inches="tight")
        return [im]

    ani = animation.FuncAnimation(fig, update, frames=grid_size)

    if slider:
        #Sliders
        plt.subplots_adjust(bottom=0.2)
        axwave = plt.axes([0.25, 0.05, 0.5, 0.03])

        sliderwave = Slider(axwave, 'Z axis', 0, 31, valinit=0, valfmt='%d')
        sliderwave.on_changed(update)

    plt.show()

#function to plot motion of particles
#uses a list of particle positions across all timesteps generated using the append_new_array function below
def plot_particle_motion(particle_positions, num_steps, interval,save_gif=False,name="particles.gif",particle_size=0.1):
    """
    Plots and shows the motion of particles over time.

    :param particle_positions:
                List of arrays containing particle positions at each timestep
                It is a (total_time_steps,num_particles,3) list
                Can be constructed using the append_new_array function
                
    :param num_steps:
                Total number of timesteps in the simulation.
    :param interval:
                Interval (in milliseconds) between frames in the animation.

    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Set up the plot elements
    ax.set_xlim3d(0, 32)
    ax.set_ylim3d(0, 32)
    ax.set_zlim3d(0, 32)
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_zlabel('Z-axis')

    #making the animated plot
    def update_graph(num):
        data = particle_positions[num]
        ax.clear()
        ax.set_xlim3d(0, 32)
        ax.set_ylim3d(0, 32)
        ax.set_zlim3d(0, 32)
        ax.scatter3D(data[..., 0], data[..., 1], data[..., 2], s=particle_size)
        ax.set_title("Particle Motion at t = {}".format(num))
        return ax

    ani = animation.FuncAnimation(fig, update_graph, num_steps,
                                             interval=interval, blit=False)

    if save_gif:
        # To save the animation using Pillow as a gif
        writer = animation.PillowWriter(fps=15,
                                        metadata=dict(artist='Me'),
                                        bitrate=1800)
        ani.save(name, writer=writer)
    else:
        plt.show()


def plot_motion_from_save(name,num_particles,save_gif=False,fname="particles.gif",particle_size=0.1):
    # Read the saved positions
    loaded_moves = np.loadtxt(name)

    # Calculate the correct number of timesteps
    num_dimensions = 3  # x, y, z coordinates
    total_elements = loaded_moves.size
    num_timesteps = total_elements // (num_particles * num_dimensions)

    # Reshape back to original shape
    loaded_moves = loaded_moves.reshape((num_timesteps, num_particles, num_dimensions))


    plot_particle_motion(loaded_moves, num_timesteps, interval = 3.2,save_gif=save_gif,name=fname,particle_size=particle_size)


# function to append a new (num_particles, 3) array to the existing 3D array
def append_new_array(existing_array, new_array):
    # reshape new array to (1, num_particles, 3)
    new_array_reshaped = new_array[np.newaxis, ...]
    return np.concatenate((existing_array, new_array_reshaped), axis=0)

# function for plotting a single slice of a density field, i.e., a 2D slice of field at a position along an axis
def plot_density_slice(density_field, time, axis='z', slice_position=16, grid_size=32, save_plot=False):
    # initialize plot
    fig, ax = plt.subplots()

    # create plot based on indicated axis and position
    if axis == 'z':
        im = ax.imshow(density_field[:, :, slice_position], cmap='viridis',norm=colors.LogNorm(vmin=0.001,vmax=np.max(density_field)))
        ax.invert_yaxis()

        #colorbar/density in log scale to better visualize the spread of data
        fig.colorbar(im,ax=ax,label='log(Density)')
        ax.set_title("Density Field Slice (Z = " + str(slice_position) + ") at t = " + str(time))
        #setting axes ranges
        ax.set_xlim(0, grid_size)
        ax.set_ylim(0, grid_size)
        #setting axes labels
        ax.set_xlabel("X axis")
        ax.set_ylabel("Y axis")
    elif axis == 'y':
        im = ax.imshow(density_field[:, slice_position, :], cmap='viridis',norm=colors.LogNorm(vmin=0.001,vmax=np.max(density_field)))
        ax.invert_yaxis()

        #colorbar/density in log scale to better visualize the spread of data
        fig.colorbar(im,ax=ax,label='log(Density)')
        ax.set_title("Density Field Slice (Y = " + str(slice_position) + ") at t = " + str(time))
        #setting axes ranges
        ax.set_xlim(0, grid_size)
        ax.set_ylim(0, grid_size)
        #setting axes labels
        ax.set_xlabel("X axis")
        ax.set_ylabel("Z axis")
    else:
        im = ax.imshow(density_field[slice_position, :, :], cmap='viridis',norm=colors.LogNorm(vmin=0.001,vmax=np.max(density_field)))
        ax.invert_yaxis()

        #colorbar/density in log scale to better visualize the spread of data
        fig.colorbar(im,ax=ax,label='log(Density)')
        ax.set_title("Density Field Slice (X = " + str(slice_position) + ") at t = " + str(time))
        #setting axes ranges
        ax.set_xlim(0, grid_size)
        ax.set_ylim(0, grid_size)
        #setting axes labels
        ax.set_xlabel("Y axis")
        ax.set_ylabel("Z axis")
    
    # save or show
    if save_plot:
        plt.savefig('density_slice_' + axis + '_equals_' + str(slice_position))
    else:
        plt.show()

# function to plot the density field at various times
# takes array of positions at times, array of desired times to plot at, desired axis, desired slice positon
# boolean indicated whether to save plots
def plot_density_at_times(positions_at_times, time_steps, axis='z', slice_position=16, save=False):
    for t in time_steps:
        density = density = cic_density(positions_at_times[t],grid_size=32)
        plot_density_slice(density, time=t, axis=axis, slice_position=slice_position, save_plot=save)
