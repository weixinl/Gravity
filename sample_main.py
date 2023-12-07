#Gravity on a mesh
import matplotlib.pyplot as plt
import numpy as np

from densities import *
from distribute import *
from plotter import *
from integrator import *
import dynamicsteps


def load_txt(name,num_particles=32**3):
    # Read the saved positions
    loaded_moves = np.loadtxt(name)

    # Calculate the correct number of timesteps
    num_dimensions = 3  # x, y, z coordinates
    total_elements = loaded_moves.size
    num_timesteps = total_elements // (num_particles * num_dimensions)

    # Reshape back to original shape
    loaded_moves = loaded_moves.reshape((num_timesteps, num_particles, num_dimensions))

    return loaded_moves,num_timesteps

def main():
    # Generate particles and density field
    grid_size = 32
    num_particles = 10

    ####---------load file---------------########
    fname = 'test.txt'
    moves, current_tstep = load_txt(fname, num_particles)
    ####------------------------------------########

    # if loaded file, comment out the distribution function stuff below

    ##Gaussian distribution parameters
    center = (grid_size/2, grid_size/2,grid_size/2)
    a = 5
    b_to_a = 0.5
    c_to_a = 0.8

    #particles = distribute_spherical(radius=10,num_particles=num_particles)
    particles = distribute_gaussian(center, a, b_to_a, c_to_a,grid_size=grid_size,num_particles=num_particles)
    #particles = np.zeros(np.shape(particles))
    # particles[0] = [14,16,16]
    # particles[1] = [18,16,16]
    # particles[2] = [16,16,20]
    # particles[3] = [16,15,15]
    # for i in range(100):
    #     particles[i] = [16,15,15]
    # particles[-1] = [14,16,16]

    density = cic_density(particles,grid_size=grid_size)
    #plot_particles(particles)
    #plot_density(density,grid_size)

    #array to store all positions at each timestep
    moves = np.array(particles)
    #adding a new axis to store timesteps: moves[0] gives positions of all particles at t=0
    moves = moves[np.newaxis,...]

    #testing out with zeroes first
    velocity = np.zeros(np.shape(particles))

    num_steps = 50
    time_step = 1

    for i in range(num_steps):
        #time_step = dynamicsteps.step_by_dist(particles)
        particles,velocity= integrate(particles,velocity,density,time_step=time_step)
        density = cic_density(particles,grid_size=grid_size)
        # enforcing toroidal boundary conditions
        particles = particles%32
        moves = append_new_array(moves,particles)

    plot_particle_motion(moves, num_steps, time_step, interval = 3.2,save_gif=False,name="32_3particles_50s.gif",particle_size=10)
    #plot_motion_from_save("test.txt",num_particles)

    #============================save file ========================================#
    #save file name
    fname = 'test.txt'

    # reshaping the array from 3D to 2D to save as txt file
    # automatically converts to appropriate shape when reading it back in later
    moves_reshaped = moves.reshape(moves.shape[0], -1)
    np.savetxt(fname, moves_reshaped)


if __name__ == '__main__':
    main()

