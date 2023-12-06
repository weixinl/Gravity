import numpy as np
import matplotlib.pyplot as plt
from distribute import *

def decreasing_step(t):
    '''
    decreasing funciton of time t
    t: current time
    return: current dynamic time step
    '''
    step_coeff = 0.1
    return step_coeff / (1 + np.exp(t))

def step_by_dist(particles):
    '''
    dynamic steps based on average distance
    small step for small distance
    particles: N * 3 array
    return: current dynamic time step
    '''
    N = particles.shape[0]
    dist = 0
    for i in range(N - 1):
        for j in range(i + 1, N):
            diff = particles[i] - particles[j]
            dist += np.sqrt(np.sum(diff * diff))
    dist_avg = 2 * dist/(N * (N - 1))
    step_coeff = 1e-2 
    return step_coeff * dist_avg    

if __name__ == "__main__":
    # steps = 100
    # step_sizes = []
    # for i in range(steps):
    #     step_sizes.append(decreasing_step(i))
    # plt.plot(step_sizes)
    # plt.show()
        # Generate particles and density field
    grid_size = 32
    num_particles = 101
    center = (grid_size/2, grid_size/2,grid_size/2)
    a = 5
    b_to_a = 0.5
    c_to_a = 0.8

    #particles = distribute_spherical(radius=10,num_particles=num_particles)
    particles = distribute_particles(center, a, b_to_a, c_to_a,grid_size=grid_size,num_particles=num_particles)
    step_size = step_by_dist(particles)
    print(step_size)