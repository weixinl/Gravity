#Gravity on a mesh
import matplotlib.pyplot as plt
import numpy as np

from densities import *
from distribute import *
from plotter import *
from integrator import *
from solver import *

def main():
    # Generate particles and density field
    grid_size = 32
    num_particles = 100
    center = (grid_size/2, grid_size/2,grid_size/2)
    a = 4
    b_to_a = 1
    c_to_a = 1

    #particles = distribute_spherical(radius=10,num_particles=num_particles)
    particles = distribute_particles(center, a, b_to_a, c_to_a,grid_size=grid_size,num_particles=num_particles)

    density = cic_density(particles,grid_size=grid_size)

    #array to store all positions at each timestep
    moves = np.array(particles)
    #adding a new axis to store timesteps: moves[0] gives positions of all particles at t=0
    moves = moves[np.newaxis,...]

    #testing out with zeroes first
    velocity = np.zeros(np.shape(particles))
    # array to store all velocities at each time step
    velocities = velocity
    velocities = velocities[np.newaxis, ...]
    # array to store potentials at each time step
    potentials = solve_poisson(density)
    potentials = potentials[np.newaxis, ...]

    tsteps = 100

    for i in range(tsteps):
        particles,velocity,potential= integrate(particles,velocity,density,time_step=0.1)
        density = cic_density(particles,grid_size=grid_size)
        # enforcing toroidal boundary conditions
        particles = particles%32
        moves = append_new_array(moves,particles)
        velocities = append_new_array(velocities, velocity)
        potentials = append_new_array(potentials, potential)


    # test virial theorem
    virial_plots(moves, velocities, potentials)



if __name__ == '__main__':
    main()


