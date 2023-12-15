#Gravity on a mesh
#plot helper

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import matplotlib.colors as colors
from matplotlib.widgets import Slider
from densities import *
import mytimer

#function to plot particles with given positions
#used for plotting initial distribution or a positions during a timestep
def plot_particles(positions,grid_size=32,save_fig=False,particle_size=0.1, title = "Distributing particles with a Gaussian Distribution", \
                   filename = "initial_distribution.png", dotsize = 0.01):
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
    plt.title(title)

    #making the axis labels closer to the axis
    ax.xaxis.labelpad = -10
    ax.yaxis.labelpad = -10
    ax.zaxis.labelpad = -10

    #making the axis ticks closer to the axis
    ax.tick_params(axis='x', pad=-5)
    ax.tick_params(axis='y', pad=-5)
    ax.tick_params(axis='z', pad=-4)

    ax.scatter3D(positions[..., 0], positions[..., 1], 0 , c='grey', s= dotsize)
    ax.scatter3D(0, positions[..., 1], positions[..., 2], c='grey', s=dotsize)
    ax.scatter3D(positions[..., 0], grid_size , positions[..., 2], c='grey', s=dotsize)

    if save_fig:
        plt.savefig(filename, dpi=1200, bbox_inches="tight")
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
def plot_particle_motion(particle_positions, timesteps, interval,save_gif=False,name="particles.gif",particle_size=0.1):
    """
    Plots and shows the motion of particles over time.

    :param particle_positions:
                List of arrays containing particle positions at each timestep
                It is a (total_time_steps,num_particles,3) list
                Can be constructed using the append_new_array function
                
    :param timer: mytimer.Timer
                object timer manages current time, dynamic time steps and time log

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

    ani = animation.FuncAnimation(fig, update_graph, timesteps,
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
        plt.savefig('sphere_density_slice_' + axis + '_equals_' + str(slice_position) + '_t_' + str(time))
    else:
        plt.show()

# function to plot the density field at various times
# takes array of positions at times, array of desired times to plot at, desired axis, desired slice positon
# boolean indicated whether to save plots
def plot_density_at_times(positions_at_times, time_steps, axis='z', slice_position=16, save=False):
    for t in time_steps:
        density = cic_density(positions_at_times[t],grid_size=32)
        plot_density_slice(density, time=t, axis=axis, slice_position=slice_position, save_plot=save)

# function for plotting potential vs. kinetic energy
# takes arrays of particle positions, particle velocities, and gravitational potentials at each time step
# takes time index to start the plot at in case you want to only plot later times
# returns arrays of U and K at each time step
def virial_plots(positions_at_times, velocities_at_times, potentials_at_times, start_time=0,save=False):
    # calculate density at each time step
    densities = []
    for pos in positions_at_times:
        densities.append(cic_density(pos,grid_size=32))
    densities = np.array(densities)

    # calculate potential energy at each position at each time step (particle mass = 1)
    # factor of 0.5 added after a discussion with Mike
    U_positions = 0.5 * densities * potentials_at_times
    # array of total potential energy at each time step
    U = []
    for u in U_positions:
        U.append(np.sum(u))
    U = np.array(U)

    # calculate kinetic energy at each time step (m=1)
    # convert velocities to speeds
    speeds_at_times = velocity_to_speed(velocities_at_times)
    # array of kinetic energies for each particle at each time
    K_particles = 0.5 * speeds_at_times**2
    # array of total kinetic energy at each time step
    K = []
    for K_time in K_particles:
        K.append(np.sum(K_time))
    K = np.array(K)
    
    fit, fit_coeffs = fit_data(K[start_time:], U[start_time:], 1)

    plt.plot(K[start_time:], U[start_time:], '.', label='Data')
    #plt.plot(K[start_time:], fit, 'k-', label='U = ' + str(fit_coeffs[0]) + 'K + ' + str(fit_coeffs[1]))
    plt.plot(K, -2*K, 'r-', label='Virial Theorem')
    plt.title('Total Potential Energy vs Total Kinetic Energy')
    plt.xlabel('Kinetic Energy')
    plt.ylabel('Potential Energy')
    plt.legend()
    if save==True:
        plt.savefig('gravity_virial_test.png')
    else:
        plt.show()
    return K, U

# function to take magnitude of velocity for each particle
# returns array of speeds for each particle at each time
def velocity_to_speed(velocities_at_times):
    speeds = []
    for vel in velocities_at_times:
        speeds.append(np.sqrt(np.sum(vel**2)))
    return np.array(speeds)

# function for creating linear fits
# create fit line
def fit_data(x, y, degree):
    fit_coeffs = np.polyfit(x, y, deg=degree)
    return fit_coeffs[0]*x + fit_coeffs[1], fit_coeffs

def virial_from_load(K, U, start_time=0, save=False):
    plt.plot(K[start_time:], U[start_time:], '.', label='Data')
    #plt.plot(K[start_time:], fit, 'k-', label='U = ' + str(fit_coeffs[0]) + 'K + ' + str(fit_coeffs[1]))
    plt.plot(K, -2*K, 'r-', label='Virial Theorem')
    plt.title('Total Potential Energy vs Total Kinetic Energy')
    plt.xlabel('Kinetic Energy')
    plt.ylabel('Potential Energy')
    plt.xlim(0.4e8, 0.75e8)
    plt.ylim(-1.5e8, -0.5e8)
    plt.legend()
    if save==True:
        plt.savefig('gravity_virial_test_zoom.png')
    else:
        plt.show()
    

    


