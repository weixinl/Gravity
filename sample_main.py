#Gravity on a mesh
from densities import *
from distribute import *
from plotter import *
import solver
from integrator import *

def main():
    # Generate particles and density field
    grid_size = 32
    num_particles = 2
    center = (grid_size/2, grid_size/2,grid_size/2)
    a = 5
    b_to_a = 0.5
    c_to_a = 0.8

    particles = distribute_particles(center, a, b_to_a, c_to_a,grid_size=grid_size,num_particles=num_particles)
    particles = np.zeros(np.shape(particles))
    particles[0] = [14,14,16]
    particles[1] = [18,18,16]
    # particles[2] = [16,16,15]
    # particles[3] = [16,15,15]

    density = cic_density(particles,grid_size=grid_size)
    # plot_particles(particles)
    # plot_density(density,grid_size)

    #array to store all positions at each timestep
    moves = np.array(particles)
    #adding a new axis to store timesteps: moves[0] gives positions of all particles at t=0
    moves = moves[np.newaxis,...]

    #testing out with random values for now
    velocity = np.zeros(np.shape(particles))

    tsteps = 100

    for i in range(tsteps):
        particles,velocity= integrate(particles,velocity,density,1)
        density = cic_density(particles,grid_size=grid_size)
        # enforcing toroidal boundary conditions
        particles = particles%32
        moves = append_new_array(moves,particles)

    # np.savetxt('particles.txt', moves.reshape(-1, moves.shape[-1]))
    plot_particle_motion(moves, tsteps, interval = 3.2,save_gif=False,particle_size=5)




if __name__ == '__main__':
    main()

