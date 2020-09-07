
def CalculatePvalues(dists, exp_inp):
    """calculating  pvalues for known level of expression/ making bh correction 
       :dists: - matrix of counts
       :exp_inp: - vector with expression level
       :return: - [Pvals_low,Pvals_up,Pvals_inter_low,Pvals_inter_up,Pvals_low_adj,Pvals_up_adj,Pvals_inter_low_adj,Pvals_inter_up_adj] 
                - pvalues
    """
    r, _= exp_inp.shape # finding the shapes of exp_inp for inizialazing arrays for pvals
    c = exp_inp[0].shape[0] # finding the shapes of exp_inp for inizialazing arrays for pvals
    #print('r',r)
    #print('c',c)
    
    Pvals_low = np.zeros((r,c))#inizialazing arrays for pvals
    Pvals_up = np.zeros((r,c))#inizialazing arrays for pvals
    Pvals_inter_low = np.zeros((r,c))#inizialazing arrays for pvals
    Pvals_inter_up = np.zeros((r,c))#inizialazing arrays for pvals
    #calculating pvals with several jobs
    Parallel(n_jobs=-2)(delayed(as_par_for)(i, exp_inp, dists,c, Pvals_low, Pvals_up, Pvals_inter_low,Pvals_inter_up) for i in range(r)) 
    #print('here -------------------> ' ,Pvals_low)
    
    #making bh correction. reshape is needed for multipletests 
    tmp = Pvals_low.reshape((r*c))
    #print('---------->', tmp.shape)
    adjp = multipletests(tmp, method = 'fdr_bh')[1]
    Pvals_low_adj = adjp[0:c*r].reshape((c,r)).T
    
    tmp = Pvals_up.reshape((c*r))
    adjp =  multipletests(tmp, method = 'fdr_bh')[1]
    Pvals_up_adj = adjp[0:c*r].reshape((c,r)).T
    
    tmp = Pvals_inter_low.reshape((c*r))
    adjp = multipletests(tmp, method = 'fdr_bh')[1]
    Pvals_inter_low_adj = adjp[0:c*r].reshape((c,r)).T
    
    tmp = Pvals_inter_up.reshape((c*r))
    adjp = multipletests(tmp, method = 'fdr_bh')[1]
    Pvals_inter_up_adj = adjp[0:c*r].reshape((c,r)).T
    
    return [Pvals_low,Pvals_up,Pvals_inter_low,Pvals_inter_up,Pvals_low_adj,Pvals_up_adj,Pvals_inter_low_adj,Pvals_inter_up_adj] 
        
        
        
            

