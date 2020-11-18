def CalculatePvalues(dists, exp_inp):
  """calculating  pvalues for known level of expression/ making bh correction 
       :dists: - matrix of counts
       :exp_inp: - vector with expression level
       :return: - [Pvals_low,Pvals_up,Pvals_inter_low,Pvals_inter_up,Pvals_low_adj,Pvals_up_adj,Pvals_inter_low_adj,Pvals_inter_up_adj] 
                - pvalues
    """
  from rewrite_MatLab.BooleanizeFromFile import as_par_for as as_par_for
  import numpy as np
  from joblib import Parallel as Parallel
  from joblib import delayed as delayed
  r, _= exp_inp.shape  #finding the shapes of exp_inp for inizialazing arrays for pvals
    #print('c',c)
  c = exp_inp[0].shape[0] # finding the shapes of exp_inp for inizialazing arrays for pvals
  print('r',r)
  xq = np.arange(0, max(exp_inp)+2*1e-3, 1e-3)
  Pvals_low = np.zeros((r,c))#inizialazing arrays for pvals
  Pvals_up = np.zeros((r,c))#inizialazing arrays for pvals
  Pvals_inter_low = np.zeros((r,c))#inizialazing arrays for pvals
  Pvals_inter_up = np.zeros((r,c))#inizialazing arrays for pvals

    #calculating pvals with several jobs
  Parallel(n_jobs=4,prefer="threads")(delayed(as_par_for)(i, exp_inp, xq, dists,c, Pvals_low, Pvals_up, Pvals_inter_low,Pvals_inter_up) for i in range(r))
  print('here -------------------> ' ,Pvals_low)

  where_are_NaNs = isnan(Pvals_low)
  Pvals_low[where_are_NaNs] = 1

    #making bh correction. reshape is needed for multipletests 
  tmp = Pvals_low.reshape((r*c))
    #print('dists_inp[137]---------->', dists[137][0].shape)
  print('is printed before adjp  ------>', Pvals_low[145:160])
  adjp = multipletests(tmp, method = 'fdr_bh')[1]
  Pvals_low_adj = adjp[0:c*r].reshape((c,r)).T#почему такой срез??
    #print(Pvals_low_adj[1:11])
    
  tmp = Pvals_up.reshape((c*r))
  adjp =  multipletests(tmp, method = 'fdr_bh')[1]
  Pvals_up_adj = adjp[0:c*r].reshape((c,r)).T
  tmp = Pvals_inter_low.reshape((c*r))
  adjp = multipletests(tmp, method = 'fdr_bh')[1]
  Pvals_inter_low_adj = adjp[0:c*r].reshape((c,r)).T
  tmp = Pvals_inter_up.reshape((c*r))
  adjp = multipletests(tmp, method = 'fdr_bh')[1]
  Pvals_inter_up_adj = adjp[0:c*r].reshape((c,r)).T
  print('is printed before exit  ------>', Pvals_low_adj[130:145])
  return [Pvals_low,Pvals_up,Pvals_inter_low,Pvals_inter_up,Pvals_low_adj,Pvals_up_adj,Pvals_inter_low_adj,Pvals_inter_up_adj] 
        
        

  
