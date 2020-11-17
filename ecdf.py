def ecdf(x, freq):
    """calculete empirical destribtion function for values x with frequences freq
       :x:list - values 
       :freq:list - frequences of x
       :return: x - x values for ecdf
       :return: ys - y values for ecdf 
    """
    import numpy as np
    #print(len(np.nonzero(x)[0]))
    if len(np.nonzero(x)[0])==0:
        lenth = len(x)
        x = np.arange(0, lenth, 1)
        ys = np.ones(lenth)
        #print('here!!!')
        return x, ys
    #print('first x ---------->', x)
    indices = x.argsort()

    x = x[indices]  # -> [1 2 3]
    freq = freq[indices]
    
    
  
    if freq[0]==0:
        
        ys = np.zeros(len(x))
        ys[1:] = cumsum(freq)/sum(freq)
    else:
        ys = np.cumsum(freq)/np.sum(freq)
    
    #print('second x --------------------->', x)
    if not isinstance(x[0], float):
        print('not float',x[0])
        x = [x[i][0] for i in range(len(x))]
    return x, ys
