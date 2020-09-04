
# coding: utf-8

# In[ ]:


def BooleanizeFromFile(InputData, OutputFile ):
    import scipy.io
    import numpy as np 
    mat = scipy.io.loadmat('Thresholds_CoTFcRF.mat')#возвращает словарь. под ключем dists "матрица": [0:2136][0][0][0:1][0:1000]
    #dt_1 = np.dtype([('Ens', 'S16'), ('Triv', 'S10')], align=True)
    #with open('background_genes.txt') as inp:
    #    backgroundgenes, triv = np.loadtxt(inp, dtype=dt_1, unpack = True)
    #print(type(backgroundgenes[1]))
    dt_1 = np.dtype([('Ens', (str, 16)), ('Triv', (str, 10))], align=True)
    with open('background_genes.txt') as inp:
        backgroundgenes = np.loadtxt(inp, dtype=dt_1)#считываем данные в numpy массив пар значений(возможно, придется исправить)
    print(InputData, OutputFile)
    dt_2 = np.dtype([('Ens', (str, 20)), ('Triv', (str,10)), ('exp', int)])#то же самое, но тройки значений
    with open(InputData) as InpD:
        inpData = np.genfromtxt(InpD, dtype = dt_2, names = True)
    
    inputmat = -1*np.ones((2137,1)) #может быть исправить длину 2137??
    #print(inpData[0])
    for i in range(2137):# тут то же самое 
        #print(backgroundgenes[i][0])
        #print(inpData['ensembl_gene_id'][1])
        idx = np.nonzero(backgroundgenes[i][0] == inpData['ensembl_gene_id'])#может быть поблемма из-за регистра
        if idx[0]:
            inputmat[i] = inpData['TPM'][idx[0][0]]
    existing = np.nonzero(inputmat>-1)[0]
    #print(existing)
    dists_inp = np.array([mat['dists'][i][0][0] for i in existing])#выглядит ужасно, надо исправить 
    backgroundgenes_inp = [backgroundgenes['Triv'][i] for i in existing]
    #print(dists_inp)
    exp_inp = [inputmat[i] for i in existing]
    exp_inp = np.array(exp_inp)
    Pvals_low,Pvals_up,Pvals_inter_low,Pvals_inter_up,Pvals_low_adj,Pvals_up_adj, Pvals_inter_low_adj,Pvals_inter_up_adj = CalculatePvalues(dists_inp, exp_inp)
    s = len(exp_inp)
    
    
    #тут я остановилась
    Bool_exp = np.zeros((s,1))
        
    
    for i in np.nonzero(Pvals_low_adj>0.1)[0]:
        Bool_exp[i]=1
    output = [(backgroundgenes_inp[gene], Bool_exp[gene]) for gene in range(s)]

