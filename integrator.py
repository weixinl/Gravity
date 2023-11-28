# Graduate Computational Physics
# Gravity on a Mesh Project
# Integrator file

import numpy as np
import densities
import solver

# function to run one time step of the integrator
def integrate(positions, velocities, potential, time_step):
    '''
    Integrates the particles' positions and velocities forward by one time step
    using the gradient of the gravitational potential and the Verlet integration method

    :param positions: array size (num_particles, 3)
                    array of position values (x, y, z) for each particle
    :param velocities: array size (num_particles, 3)
                    array of velocity values (vx, vy, vz) for each particle
    :param potential: array size (32, 32, 32)
                    array of potential values at each point in 32x32x32 space
    :param time_step: float
                    size of whole time step
    :return: new_positions:
                    returns array of positions values (x, y, z) after one time step
    :return: new_velocities:
                    returns array of positions values (x, y, z) after one time step
    :return: new_potential:
                    returns array of potential values
    '''
    
    # get values of force (3 components) on each particle at the current positions
    F_current = grav_force(potential, positions)

    # calculate v(t + half step)
    v_half = velocities + 0.5*time_step*F_current
    # calculate x(t + step)
    new_positions = positions + time_step*v_half

    # calculate updated density with new positions
    new_density = densities.cic_density(new_positions)
    # update potential
    new_potential = solver.solve_poisson(new_density)

    # calculate new values of force from the new potential
    new_F = grav_force(new_potential, new_positions)
    # calculate v(t + step)
    new_velocities = v_half + 0.5*time_step*new_F

    # return position, velocity, new_potential after a whole time step
    return new_positions, new_velocities, new_potential

# function that calculates the gravitational force on each particle
def grav_force(potential, positions):
    '''
    Calculates the gravitational force on each particle
    Takes the gradient of the potential and returns the value of the potential gradient (force) nearest
    to this position (approximation since our potential is a discrete function)

    :param potential: array size (32, 32, 32)
                    array of potential values at each point in 32x32x32 space
    :param positions: array size (num_particles, 3)
                    array of position values for each particle (x, y, z)
    :return: F: array size (num_particles, 3)
                    returns array of force on each particle in each direction (Fx, Fy, Fz)
    '''
    # take gradient of the potential
    dphi = np.gradient(potential)

    # now we have values of force at discrete position values on the grid
    # now I am going to round every particle position to the nearest whole number
    # in order to get the closest approximation of the force on each particle
    rounded_positions = np.clip(np.rint(positions), 0, 31)
    # convert these to integers to be used as indices of the potential array
    int_positions = rounded_positions.astype(int)

    # create arrays of position values in each dimension
    xpos = int_positions[:,0]
    ypos = int_positions[:,1]
    zpos = int_positions[:,2]

    # make an array of force values at each particle's position
    # F is array of size (num_particles, 3) with force values (Fx, Fy, Fz) for each particle at their
    # current position
    F = [[0]*3]*np.size(xpos)
    for i in range(np.size(xpos)):
        # Fx = dphi/dx
        F[i][0] = dphi[0][xpos[i], ypos[i], zpos[i]]
        # Fy = dphi/dy
        F[i][1] = dphi[1][xpos[i], ypos[i], zpos[i]]
        # Fz = dphi/dz
        F[i][2] = dphi[2][xpos[i], ypos[i], zpos[i]]
        
    
    # return F as an array
    return np.array(F)

# placeholder for laplace equation solver (will come from different source file)
#def solve_laplace(density):
#    potential = density
#    return potential

