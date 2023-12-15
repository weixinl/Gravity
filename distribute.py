#Gravity on a mesh
#distribute particles with multivariate Gaussian

import numpy as np
import utils

def distribute_gaussian(center=(16.,16.,16.), a=5., b_to_a=0.5, c_to_a=0.8, grid_size =32, num_particles=32 ** 3):
    '''
    Distribute particles according to a multivariate Gaussian.
    :param center: tuple
                    center of gaussian
    :param a: np.float
                    semi-major axis
    :param b_to_a: np.float
                    ratio b/a (b is  semi minor axis along y)
    :param c_to_a: np.float
                    ratio c/a (c is semi minor axis along z)
    :param grid_size: int
                    length of each axis of the grid
    :param num_particles: int
                    number of particles (deafult is 32^3)
    :return: positions: np.array
                    returns array of shape (num_particles,3) containing the x,y,z coordinates of each particle
    '''

    mean = np.array(center)

    # Covariance matrix, aligned with the coordinate axes.
    #https://cs229.stanford.edu/section/gaussians.pdf
    cov = np.diag([a**2 , (b_to_a * a)**2, (c_to_a * a)**2 ])
    
    # Generate multivariate Gaussian distributed points.
    positions = np.random.multivariate_normal(mean, cov, num_particles)

    # enforcing toroidal boundary conditions
    positions = positions%32

    return positions

def distribute_spherical(center=(16,16,16),radius=10,grid_size =32, num_particles=32 ** 3):
    #r values
    u = np.random.uniform(0.0, 1, (num_particles))
    R = radius * (u)**(1/3)

    #phi values
    phi = np.random.uniform(0.0,2*np.pi,(num_particles))

    #theta values
    costheta = np.random.uniform(-1,1,(num_particles))
    theta = np.arccos(costheta)

    x = R * np.sin(theta) * np.cos(phi) + (center[0])
    y = R * np.sin(theta) * np.sin(phi) + (center[1])
    z = R * np.cos(theta) + (center[2])

    positions = np.stack((x,y,z),axis=1)

    # enforcing toroidal boundary conditions
    positions = positions % 32

    return positions

def distribute_circle(direction=(1, 0, 0), center=(16,16,16),radius=10,grid_size =32, num_particles=20):
    '''
    an evenly distributed circle

    :param direction: tuple
                    direction of the circle (perpendicular to the circle plane)
    :param center: tuple
                    center of the circle
    :param radius: int
                    radius of the circle
    :param grid_size: int
                    length of each axis of the grid
    :param num_particles: int
                    number of particles (deafult is 32^3)
    :return: positions: np.array
                    returns array of shape (num_particles,3) containing the x,y,z coordinates of each particle
    '''

    # convert direction into np array
    direction = np.asarray(direction)

    # normalize direction to a unit vector
    direction = utils.normalize_array(direction)

    # A plane is determined by the direction vector, just direction vector cross product another vector (vec_aux) (not parallel)
    # so rotate direction vector to get another vector (vec_aux)
    vec_aux = utils.rotate_x(direction, np.pi/2)
    if(utils.calc_array_len(direction - vec_aux) < 1e-8):
        # parallel, meaningless to cross product
        vec_aux = utils.rotate_y(direction, np.pi/2)
    
    # print("vec_aux", vec_aux)

    vec_list = np.zeros((num_particles, 3))

    # cross product to get a vector in the expected plane (perpendicular to direction vector)
    vec0 = np.cross(direction, vec_aux)

    # print("vec0", vec0)

    # adjust the length of vec0
    vec0 = utils.normalize_array(vec0) * radius


    theta = np.pi * 2 / num_particles

    vec_list[0, :] = vec0
    for i in range(1, num_particles):
        vec_list[i, :] = utils.rotate_by_axis(vec0, direction, theta * i)

    # for each vector, add center position to get the particle positions
    positions = np.zeros((num_particles, 3))
    for i in range(num_particles):
        positions[i, :] = vec_list[i, :] + center

    return vec_list, positions
    