
# coding: utf-8

# In[ ]:


def ecdf(x, freq):
    
    x, freq  = zip(*sorted(zip(x,freq)))
    #xs = np.sort(x)
    if freq[0]==0:
        ys = np.zeros(len(x))
        ys[1:] = np.cumsum(freq)[1:]/np.sum(freq)
    else:
        ys = np.cumsum(freq)/np.sum(freq)
    #ys = np.arange(1, len(xs)+1)/float(len(xs))
    return x, ys

