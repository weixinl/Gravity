#Gravity on a mesh
#distribute particles with multivariate Gaussian

import numpy as np


def distribute_particles(center=(16.,16.,16.), a=5., b_to_a=0.5, c_to_a=0.8, grid_size =32, num_particles=32 ** 3):
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
    print(cov.shape)
    print(cov)
    # Generate multivariate Gaussian distributed points.
    positions = np.random.multivariate_normal(mean, cov, num_particles)

    # Clip the particles to ensure they are within the grid boundaries
    positions = np.clip(positions, 0, grid_size)

    return positions

