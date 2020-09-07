
def ecdf(x, freq):
    """calculete empirical destribtion function for values x with frequences freq
       :x:list - values 
       :freq:list - frequences of x
       :return: x - x values for ecdf
       :return: ys - y values for ecdf 
    """
    #print(len(np.nonzero(x)[0]))
    if len(np.nonzero(x)[0])==0:# if all x are 0 
        x = np.arange(0, len(x), 1)
        ys = [1 for i in range(len(x))]
        return x, ys
    #print('first x ---------->', x)
    x, freq  = zip(*sorted(zip(x,freq)))
    #xs = np.sort(x)
    if freq[0]==0:
        ys = np.zeros(len(x))
        ys[1:] = np.cumsum(freq)[1:]/np.sum(freq)
    else:
        ys = np.cumsum(freq)/np.sum(freq)
    #ys = np.arange(1, len(xs)+1)/float(len(xs))
    #print('second x --------------------->', x)
    if not isinstance(x[0], float):
        x = [x[i][0] for i in range(len(x))]
    return x, ys
