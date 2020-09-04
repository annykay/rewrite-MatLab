
# coding: utf-8

# In[ ]:


def CalculatePvalues(dists, exp_inp):
    import numpy as np
    c, _= exp_inp.shape
    #print('c',c)
    r = exp_inp[0].shape[0]
    #print('r',r)
    Pvals_low = np.zeros((r,c))
    Pvals_up = np.zeros((r,c))
    Pvals_inter_low = np.zeros((r,c))
    Pvals_inter_up = np.zeros((r,c))
    Parallel(n_jobs=-2)(delayed(as_par_for)(i, exp_inp, dists,c, Pvals_low, Pvals_up, Pvals_inter_low,Pvals_inter_up) for i in range(r))
    print('here -------------------> ' ,Pvals_low )
    tmp = Pvals_low.reshape((r*c))
    #print('---------->', tmp.shape)
    adjp = multipletests(tmp, method = 'fdr_bh')[1]
    Pvals_low_adj = adjp[0:c*r].reshape((r, c))#почему такой срез??
    tmp = Pvals_up.reshape((c*r))
    adjp =  multipletests(tmp, method = 'fdr_bh')[1]
    Pvals_up_adj = adjp[0:c*r].reshape((r,c))
    tmp = Pvals_inter_low.reshape((c*r))
    adjp = multipletests(tmp, method = 'fdr_bh')[1]
    Pvals_inter_low_adj = adjp[0:c*r].reshape((r,c))
    tmp = Pvals_inter_up.reshape((c*r))
    adjp = multipletests(tmp, method = 'fdr_bh')[1]
    Pvals_inter_up_adj = adjp[0:c*r].reshape((r,c))
    return [Pvals_low,Pvals_up,Pvals_inter_low,Pvals_inter_up,Pvals_low_adj,Pvals_up_adj,Pvals_inter_low_adj,Pvals_inter_up_adj] 
        
        
            

