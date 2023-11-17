#Gravity on a mesh

import numpy as np


def cic_density(positions, grid_size=32,cell_size =1,mass=1. ):
    '''
    Compute the density field using Cloud-in-Cell (CIC) approach
    :param positions: np.array
            array of shape (num_particles,3) containing the x,y,z coordinates of each particle
    :param grid_size: int
            number of cells in each axis of the grid
    :param cell_size: int
            length of each cubical cell, also length of particle "cloud"
    :param mass: np.float
            mass of each particle
    :return density: np.array
            returns array of shape (grid_size,grid_size,grid_size) with density field due to all the particles
     --------------------------------------------------------------------------------------------------------
    Note: we only take weights for cells next to the primary cell (x_p) since weighting is 0 for cells where delta> cell_size
          weights[0] contains weights for lower cell (cell where the particle "cloud" begins)
          and and weights[1] contains weight for higher cell (the next cell)
    '''
    

    #empty grid to store densities for each cell
    density = np.zeros((grid_size, grid_size, grid_size))

    #populating density grid according to CIC interpolation
    for x in positions:
        #indices of the lower bound of the cell containing the particle (x_p,y_p,z_p)
        x_p = np.int_(np.floor(x))

        #distance of the particle from the lower corner (|x-x_p|,|y-y_p|.|z-z_p|)
        delta = np.abs(x - x_p)

        #weights for the lower bound of the cell and all surrounding cells
        weights = [cell_size - delta, delta]

        # Add the weighted density to the lower bound and all surrounding cells
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    ii, jj, kk = (x_p[0] + i), (x_p[1] + j), (
                                x_p[2] + k)
                    #if any index goes out of bound, squishes the cloud to fit to the last cell
                    if ii>31:
                        ii=31
                    if jj>31:
                        jj=31
                    if kk>31:
                        kk=31

                    density[ii, jj, kk] += mass * weights[i][0] * weights[j][1] * weights[k][2]
    return density
