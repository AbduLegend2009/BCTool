import numpy as np
def seed_generator(X:np.ndarray):
    n,m = X.shape
    selected_genes=np.random.choice(n, size=2,replace=False)
    selected_conditions=np.random.choice(m, size=2,replace=False)
    return selected_genes, selected_conditions
def purging(X:np.ndarray , selected_genes , selected_conditions , msr_threshold, a):
    sub=X.np[A[np.ix_(selected_genes, selected_conditions)]]
    selected_gene_mean=X.mean(axis=1,keepdims=True)
    selected_conditions_mean=X.mean(axis=0,keepdims=True)
    sub_mean=sub.mean()
    residues=sub-selected_gene_mean-selected_conditions_mean+sub_mean
    msr=(residues**2).mean()
    if msr<=msr_threshold:
        return selected_genes, selected_conditions
    violations=np.abs(residues)>msr_threshold
    gene_violations=violations.sum(axis=1)
    condition_violations=violations.sum(axis=0)

    
   
