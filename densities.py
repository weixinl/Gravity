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
    for x_p in positions:
        #indices of the cell containing the particle (x_p,y_p,z_p)
        indices = np.int32(np.floor(x_p)) % grid_size

        i,j,k = indices[0],indices[1],indices[2]

        #indices of the next cell in each direction (i+1,j+1,k+1) with boundary conditions enforced
        ii, jj, kk = (i + 1) % grid_size, (j + 1) % grid_size, (k + 1) % grid_size

        #centers of the cells containing the particle (x_c,y_c,z_c)
        x_c = np.array([i+(cell_size/2),j+(cell_size/2),k+(cell_size/2)])

        #distance of the particle from the cell center (|x_p-x_c|,|y_p-y_c|.|z_p-z_c|)
        # effectively gives us the fraction of the particle in neighboring cells
        d = np.abs(x_p - x_c)

        # fraction of particle in the parent cells
        t = np.abs(cell_size - d)


        #implemntation of CIC interpolation from "Writing a PM code" by Andrey Kravtsov
        density[i,j,k] += mass * (t[0] * t[1] * t[2])
        density[ii,j,k] += mass * (d[0] * t[1] * t[2])
        density[i,jj,k] += mass * (t[0] * d[1] * t[2])
        density[i,j,kk] += mass * (t[0] * t[1] * d[2])
        density[ii,jj,k] += mass * (d[0] * d[1] * t[2])
        density[ii,j,kk] += mass * (d[0] * t[1] * d[2])
        density[i,jj,kk] += mass * (t[0] * d[1] * d[2])
        density[ii,jj,kk] += mass * (d[0] * d[1] * d[2])

    return density
