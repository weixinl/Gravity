import numpy as np


def declining_step(t):
    '''
    decreasing funciton of time t
    t: current time
    return: current dynamic time step
    '''
    return 0.1 / (1 + np.exp(t))

# if __name__ == "__main__":
#     '''
#     '''