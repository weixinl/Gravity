#Gravity on a mesh
import matplotlib.pyplot as plt
import numpy as np

from densities import *
from distribute import *
from plotter import *
from integrator import *
import dynamicsteps

def main():
    # Generate particles and density field
    grid_size = 32
    num_particles = 101
    center = (grid_size/2, grid_size/2,grid_size/2)
    a = 5
    b_to_a = 0.5
    c_to_a = 0.8

    #particles = distribute_spherical(radius=10,num_particles=num_particles)
    particles = distribute_particles(center, a, b_to_a, c_to_a,grid_size=grid_size,num_particles=num_particles)
    particles = np.zeros(np.shape(particles))
    # particles[0] = [14,16,16]
    # particles[1] = [18,16,16]
    # particles[2] = [16,16,20]
    # particles[3] = [16,15,15]
    for i in range(100):
        particles[i] = [16,15,15]
    particles[100] = [14,16,16]

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
        particles,velocity= integrate(particles,velocity,density,time_step=0.1)
        density = cic_density(particles,grid_size=grid_size)
        # enforcing toroidal boundary conditions
        particles = particles%32
        moves = append_new_array(moves,particles)

    plot_particle_motion(moves, tsteps, interval = 3.2,save_gif=False,name="3particles.gif",particle_size=10)
    #plot_motion_from_save("test.txt",num_particles)



if __name__ == '__main__':
    main()

