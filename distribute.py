#Gravity on a mesh
#distribute particles with multivariate Gaussian

import numpy as np

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