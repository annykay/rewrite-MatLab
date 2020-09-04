
# coding: utf-8

# In[ ]:


def as_par_for(r, exp_inp, dists,c,Pvals_low, Pvals_up, Pvals_inter_low,Pvals_inter_up):
   
    for i in range(r):
        xq = np.arange(0, max(exp_inp), 1e-3)
        if dists[i][0].size==0:
            Pvals_low[i,:] = np.nan 
            Pvals_up[i, :] = np.nan
            continue
        cols = dists[i][0].size
        print(cols)
        if cols > 1:
            [f_low, x_low] = ecdf(dists[i][0][0], dists[i][0][1])
            [f_up, x_up] = ecdf(dists[i][1][0], dists[i][1][1])
        else:
            [f_low, x_low] = ecdf(dists[i][0][0],[1 for i in range(len(dists[i][0][0]))])
            [f_up, x_up] = ecdf(dists[i][1][0], [1 for i in range(len(dists[i][1][0]))])
            
        x_low[0] = 0
        x_up[0] = 0
        
        if x_low[-1]<xq[-1]:
            x_low = np.append(x_low,xq[-1])
            f_low = np.append(f_low,1)
            
        if x_up[-1]<xq[-1]:
            x_up = np.append(x_up,xq[-1])
            f_up = np.append(f_up,1)
        
        if x_low[0]==x_low[1]:
            x_low = [x_low[1], x_low[2:-1]]
            f_low = [f_low[1], f_low[2:-1]]
            
        if x_up[0]==x_up[1]:
            x_up = [x_up[1], x_up[2:-1]]
            f_up = [f_up[1], f_up[2:-1]]
         
        func_q_low = interpolate.interp1d(x_low, f_low,fill_value="extrapolate")
        fq_low = func_q_low(xq)
        
        func_q_up = interpolate.interp1d(x_up, f_up,fill_value="extrapolate")
        fq_up = func_q_up(xq)
        
        low_pval_i = np.zeros(cols)
        up_pval_i = np.zeros(cols)
        inter_low_pval_i = np.zeros(cols)
        inter_up_pval_i = np.zeros(cols)
        for j in range(c):
           
            #inter_low_pval_i[0][j] = 1-fq_low[np.nonzero(xq<=exp_inp[i][j])][-1]
            #print('here ----->   ',inter_low_pval_i[0][j])
            inter_low_pval_i[j] = 1-[fq_low[i] for i in np.nonzero(xq<=exp_inp[i][j])][0][-1]
            #inter_up_pval_i[0][j] = fq_up[np.nonzero(xq>=exp_inp[i][j])][0]
            inter_up_pval_i[j] = 1-[fq_up[i] for i in np.nonzero(xq>=exp_inp[i][j])][0][0]
            #low_pval_i[0][j] = fq_low[np.nonzero(xq>=exp_inp[i][j])][0]
            low_pval_i[j] = [fq_low[i] for i in np.nonzero(xq>=exp_inp[i][j])][0][0]
            #up_pval_i[0][j] = 1- fq_up[np.nonzero(xq>=exp_inp[i][j])][-1]
            up_pval_i[j] = 1-[fq_up[i] for i in np.nonzero(xq<=exp_inp[i][j])][0][-1]
        
        #print('here -------------------> ' ,Pvals_low[i][:] ,'end here ------------->',low_pval_i)
        Pvals_low[i,:] = low_pval_i
        Pvals_up[i,:] = up_pval_i
        Pvals_inter_low[i,:] = inter_low_pval_i
        Pvals_inter_up[i,:] = inter_up_pval_i
        
        
            

