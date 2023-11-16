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

    density_field = cic_density(particles,grid_size=grid_size)
    plot_particles(particles)
    plot_density(density_field,grid_size)






if __name__ == '__main__':
    main()

