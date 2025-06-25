import numpy as np
def seed_generator(X:np.ndarray):
    n,m = X.shape
    selected_genes=np.random.choice(n, size=2,replace=False)
    selected_conditions=np.random.choice(m, size=2,replace=False)
    overall_mean=X.mean()
    return selected_genes, selected_conditions, overall_mean
def purging(X:np.ndarray,selected_genes ,selected_conditions, threshold, a, overall_mean):
    for i in range(len(selected_genes)):
        for j in range(len(selected_conditions)):
            selected_genes_mean=X[i].mean()
            selected_conditions_mean=X.T[j][i].mean()
   
