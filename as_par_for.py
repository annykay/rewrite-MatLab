def as_par_for(r, exp_inp, dists,c,Pvals_low, Pvals_up, Pvals_inter_low,Pvals_inter_up):
   """function for multiprocessing. Should be called from CalculatePvalues"""
   
    print('r ------->', r)
    
    xq = np.arange(0, max(max(exp_inp)+2*1e-3), 1e-3)# x-ax for future distribution function
    #print(len(xq))
    depth, width = dists[r][0].shape  # need for understandibg, if there are freqences
    #print(depth, width)
    
  
    #print('depth', depth)
    if depth==0:
        Pvals_low[r,:] = np.nan 
        Pvals_up[r, :] = np.nan
        return 
    
    #print('dists1 ----->', (dists[r][0][0],'dists2 ----->', dists[r][0][1]))
    if width>1:
        #calculete empirical destrebution function with frequences 
        x_low, f_low = ecdf([dists[r][0][i][0] for i in range(depth)], [dists[r][0][i][1] for i in range(depth)])
        x_up, f_up = ecdf([dists[r][1][i][0] for i in range(depth)], [dists[r][1][i][1] for i in range(depth)])
        
    else:
        #calculete empirical destrebution function without frequences 
        #print('this message is printed before ecdf else, low works')
        x_low, f_low = ecdf([dists[r][0][i] for i in range(depth)],[1 for i in range(depth)])
        #print('this message is printed before ecdf else, up works')
        x_up, f_up = ecdf([dists[r][1][i] for i in range(depth)],[1 for i in range(depth)])
    #print('this message is printed before x_low list works')

    
  
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
    
    if x_low[0]==x_low[1]:
        x_low = x_low[1:]
        f_low = f_low[1:]
    #print('len', len(x_low))    
    
    
    #print('this message is printed before 4th if works')   
    
    if x_up[0]==x_up[1]:
        x_up = x_up[1:]
        f_up = f_up[1:]
        
   
    
    
    #print('this message is printed before 1st interpolation works') 
    #print('len', len(x_low))
    #print('len', len(f_low))
    #print('x_low', x_low)
    #print('f_low', f_low)
      
    #linear interpolation for destribution finction 
    func_q_low = interpolate.interp1d(x_low, f_low, fill_value="extrapolate")
    fq_low = func_q_low(xq)
    
    #print('fq_low --->', fq_low)        
    #print('this message is printed before 2nd interpolation works')
    #print('len', len(x_low))
    #print('len f', len(f_low))
      
    #linear interpolation for destribution finction   
    func_q_up = interpolate.interp1d(x_up, f_up, fill_value="extrapolate")
    fq_up = func_q_up(xq)
        
    low_pval_i = np.zeros((1,c))
    up_pval_i = np.zeros((1,c))
    inter_low_pval_i = np.zeros((1,c))
    inter_up_pval_i = np.zeros((1,c))
    #print('r ------>', r)
    for j in range(c):        
        inter_low_pval_i[j] = 1-fq_low[np.nonzero(xq<=exp_inp[r][j])[0][-1]]
        inter_up_pval_i[j] = fq_up[np.nonzero(xq>=exp_inp[r][j])[0][0]]
        low_pval_i[j] = fq_low[np.nonzero(xq>=exp_inp[r][j])[0][0]] 
        up_pval_i[j] = 1-fq_up[np.nonzero(xq<=exp_inp[r][j])[0][-1]]
        
    #print('here -------------------> ' ,Pvals_low[r][:] ,'end here ------------->',low_pval_i)
    #print(len(Pvals_low[r]))
    Pvals_low[r][:] = low_pval_i
    Pvals_up[r][:] = up_pval_i
    Pvals_inter_low[r][:] = inter_low_pval_i
    Pvals_inter_up[r][:] = inter_up_pval_i
        
        
            
