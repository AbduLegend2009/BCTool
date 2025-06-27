import numpy as np
import random
def random_condition_selection(X:np.ndarray):
    m,n = X.shape
    random_condition=random.randint(0,n-1)
    return X[:,random_condition]

