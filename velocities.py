#Gravity on a mesh
#assign velocities

import numpy as np
import distribute
import plotter
import utils

def circle_no_evolution(direction=(1, 0, 0), radius=10,center=(16,16,16),num_particles = 20):
    '''
    Simplified version of Section 5 Question 4
    '''
    vec_list, positions = distribute.distribute_circle(direction= direction, radius = radius, center=center, num_particles=num_particles)
    plotter.plot_particles(positions, save_fig= True, title= "circle distribution", 
                           filename = "imgs/circle-distribution.png", particle_size = 10)
    G = 0.5
    v_val = np.sqrt(G * num_particles / radius) # magnitude of speed (cosmic velocity)
    velocities = np.zeros((num_particles, 3))
    for i in range(num_particles):
        v_direction = utils.normalize_array(np.cross(direction, vec_list[i]))
        velocities[i, :] = v_direction * v_val

    return positions, velocities

    

def sperical_no_evolution(radius=10,center=(16,16,16)):
    '''
    Section 5 Question 4
    For an initially spherical distribution, assign initial velocities that will lead to a nearly no evolution

    return particle positions and velocities
    '''

    # circles_directions = np.asarray([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    circles_directions = np.asarray(utils.fibonacci_sphere(10))
    circle_num = circles_directions.shape[0]

    # particle num per circle
    circle_particle_num = 20

    num_particles = circle_num * circle_particle_num


    G = 0.45
    v_val = np.sqrt(G * num_particles / radius) # magnitude of speed (cosmic velocity)

    positions_list = []
    velocities_list = []

    for circ_i in range(circle_num):

        circ_vec_list, circ_positions = distribute.distribute_circle\
            (direction= circles_directions[circ_i], radius = radius, center=center, num_particles=circle_particle_num)

        circ_velocities = np.zeros((circle_particle_num, 3))
        for i in range(circle_particle_num):
            v_direction = utils.normalize_array(np.cross(circles_directions[circ_i], circ_vec_list[i]))
            circ_velocities[i, :] = v_direction * v_val
        
        positions_list.append(circ_positions)
        velocities_list.append(circ_velocities)

    positions_all = np.concatenate(positions_list)
    velocities_all = np.concatenate(velocities_list)

    return positions_all, velocities_all

def angular_momentum(particles,velocity,center=(16,16,16),radius=10,fraction=1,mass=1):
    average_momentum = mass * np.average(velocity)

    #-------25 seems to work for 32**2 particles, have to check for other numbers-----#
    if average_momentum==0:
        average_momentum=25

    #unit vector in the axis of rotation
    z_hat = np.array([0, 0, 1])

    #vector from center of distribution to particle
    r = (particles - center)

    #magnitude of r
    r_mag = np.linalg.norm(r)

    #unit vector in r direction
    r_hat = r/r_mag

    #unit vector in the velocity direction
    #v_hat is perpendicular to r and z (i.e. will be tangential component of velocity, giving the sphere a spin)
    v_hat = np.cross(r_hat, z_hat)

    #magnitude of velocity
    # (fraction) of the (average_momentum) per particle times the (radius) i.e. characteristic size of the system
    v_mag = fraction * average_momentum * radius

    #total velocity in the v_hat direction
    velocity = v_mag *v_hat

    return velocity

# def no_movement(particles,velocity,center=(16,16,16),mass=1):
#
#     tot_mass = len(particles) * mass
#
#     # vector from center of distribution to particle
#     r = (particles - center)
#
#     # magnitude of r
#     r_mag = np.linalg.norm(r)
#
#     # unit vector in r direction
#     r_hat = r / r_mag
#
#
#     #velocities in the rhat direction to balance inward pull
#     velocity = np.sqrt(tot_mass/np.abs(r_mag))* r_hat
#
#
#     return velocity

if __name__ == "__main__":
    circle_no_evolution()
