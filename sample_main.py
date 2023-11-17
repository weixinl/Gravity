#Gravity on a mesh
from densities import *
from distribute import *
from plotter import *
from integrator import *

def main():
    # Generate particles and density field
    grid_size=32
    num_particles= 32**3
    center = (grid_size/2, grid_size/2,grid_size/2)
    a = 5
    b_to_a = 0.5
    c_to_a = 0.8

    particles = distribute_particles(center, a, b_to_a, c_to_a,grid_size=grid_size,num_particles=num_particles)

    density = cic_density(particles,grid_size=grid_size)
    plot_particles(particles)
    plot_density(density,grid_size)

    #array to store all positions at each timestep
    moves = np.array(particles)
    #adding a new axis to store timesteps: moves[0] gives positions of all particles at t=0
    moves = moves[np.newaxis,...]

    #testing out with random values for now
    velocity = np.zeros(np.shape(particles))
    potential = np.zeros(np.shape(density))
    tsteps=50

    for i in range(tsteps):
        particles,velocity,potential = integrate(particles,velocity,potential,0.1)
        particles = particles%32
        moves = append_new_array(moves,particles)


    plot_particle_motion(moves, tsteps, interval = 3.2,save_gif=False)
    






if __name__ == '__main__':
    main()

