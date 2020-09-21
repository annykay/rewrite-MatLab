def as_par_for(r, exp_inp, xq, dists,c,Pvals_low, Pvals_up, Pvals_inter_low,Pvals_inter_up):
  """function for multiprocessing. Should be called from CalculatePvalues"""
    print('r ------->', r)
    
    
    #print(len(xq))
    depth, width = dists[r][0].shape # need for understandibg, if there are freqences
    print(depth, width)
    
    #have_freq = np.nonzero(dists[r][0])[0].size
    #have_values = np.nonzero(dists[r][1])[0].size
    #print('have_freq----->', have_freq)
    #print('depth', depth)
    if depth==0:
        Pvals_low[r,:] = np.nan 
        Pvals_up[r, :] = np.nan
        print('rfor nan ---->', r)
        return 
    #[dists[r][i][0] for i in range(len(dists[r]))]
    #print('dists1 ----->', (dists[r][0][0],'dists2 ----->', dists[r][0][1]))
    if width>1:
      #calculete empirical destrebution function with frequences
        #print(dists_inp[r][0])
        #print(dists_inp[r])
        #print([dists[r][0][i][0] for i in range(depth)])
        
        print('dists ----->', dists[r][0][0][1])
        
        x_low, f_low = ecdf(dists[r][0][:, 0], dists[r][0][:, 1])
        x_up, f_up = ecdf(dists[r][1][:, 0], dists[r][1][:,1])
        #print('hi')
    else:
      #calculete empirical destrebution function without frequences 
        default = np.ones(depth)
        #print('this message is printed before ecdf else, low works')
        x_low, f_low = ecdf(dists[r][0][:,0],default)
        #print('this message is printed before ecdf else, up works')
        x_up, f_up = ecdf(dists[r][1][:, 0],default)
    #print('this message is printed before x_low list works')

    #x_low = [x_low[i][0] for i in range(len(x_low))]
    #x_up = [x_up[i][0] for i in range(len(x_up))]
    #x_low = list(x_low)
    #x_up = list(x_up)
    
    x_low[0] = 0
    x_up[0] = 0
    #print('this message is printed before 1st if works')    
    if x_low[-1]<xq[-1]:
        x_low = np.append(x_low,xq[-1])
        f_low = np.append(f_low,1)
    #print('len', len(x_low))
    #print('this message is printed before 2nd if works') 
    
    if x_up[-1]<xq[-1]:
        x_up = np.append(x_up,xq[-1])
        f_up = np.append(f_up,1)
    #print('len', len(x_low))
    #print('this message is printed before 3rd if works') 
    
    #if x_low[0]==x_low[1]:
        #x_low = x_low[1:]
        #f_low = f_low[1:]
    #print('len', len(x_low))    
    
    
    #print('this message is printed before 4th if works')   
    
    #if x_up[0]==x_up[1]:
        #x_up = x_up[1:]
        #f_up = f_up[1:]
        
   
    
    
   
    func_q_low = interpolate.interp1d(x_low, f_low, fill_value="extrapolate")
    fq_low = func_q_low(xq)
    
    
    func_q_up = interpolate.interp1d(x_up, f_up, fill_value="extrapolate")
    fq_up = func_q_up(xq)
        
    low_pval_i = np.zeros((1,c))
    up_pval_i = np.zeros((1,c))
    inter_low_pval_i = np.zeros((1,c))
    inter_up_pval_i = np.zeros((1,c))
    
    for j in range(c):
           
         
        inter_low_pval_i[j] = 1-fq_low[np.nonzero(xq<=exp_inp[r][j])[0][-1]]
        inter_up_pval_i[j] = fq_up[np.nonzero(xq>=exp_inp[r][j])[0][0]]
        low_pval_i[j] = fq_low[np.nonzero(xq>=exp_inp[r][j])[0][0]]         
        up_pval_i[j] = 1-fq_up[np.nonzero(xq<=exp_inp[r][j])[0][-1]]
   
    Pvals_low[r][:] = low_pval_i
    Pvals_up[r][:] = up_pval_i
    Pvals_inter_low[r][:] = inter_low_pval_i
    Pvals_inter_up[r][:] = inter_up_pval_i
    print('is printed in the cycle  ------>', Pvals_low[r])
        
        
            
