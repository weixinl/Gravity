import numpy as np
import solver

# test
def test():
    a = np.zeros((2,2,2))
    a[0,0,0] = 1
    print(solver.solve_poisson(a))