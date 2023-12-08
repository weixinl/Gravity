#Gravity on a mesh
#assign velocities

import numpy as np

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