def BooleanizeFromFile(InputData, OutputFile ):
   """ function for 
        :InputData: str - the name of file with input data  
        :OutputFile: str - the name of a file, where the results would be written 
        :return: -
        """
    import scipy.io#need for loading .mat files 
    import numpy as np #need for calculations 
    mat = scipy.io.loadmat('Thresholds_CoTFcRF.mat')  # load data from 'Treshholds_CoTFcRF.mat'. Returns a dict, numeric data is under 
                                                     # the key 'dists'. mat['dists'] is a tensor with the shapes : [0:2136][0][0][0:1][0:999]
    
    dt_1 = np.dtype([('Ens', (str, 16)), ('Triv', (str, 10))], align=True)# special datatype for loading data from 'background_genes.txt' 
    with open('background_genes.txt') as inp:
        backgroundgenes = np.loadtxt(inp, dtype=dt_1)#считываем данные в numpy массив пар значений(возможно, придется исправить)
    print(InputData, OutputFile)
    dt_2 = np.dtype([('ensembl_gene_id', (str, 20)), ('hgnc_symbol', (str,10)), ('TPM', int)])#то же самое, но тройки значений
    with open(InputData) as InpD:#special datatype for loading data from InputData
        inpData = np.genfromtxt(InpD, dtype = dt_2, names = True)
    
    inputmat = -1*np.ones((2137,1)) #making an array to get rid of genes, that are not in backgroundgenes
    #maybe, it's worth to write len(inpData) instead of 2137
    #print(inpData[0])
    for i in range(2137): #getting rid of genes, that are not in backgroundgenes array
        #print(backgroundgenes[i][0])
        #print(inpData['ensembl_gene_id'][1])
        idx = np.nonzero(inpData['ensembl_gene_id']==backgroundgenes['Ens'][i] )#findibg indexis of genes  in inpData, that are not in backgroundgenes

        #print(idx[0])
        if len(idx[0])>0: # here and after: np.nonzero returns two arrays: colmn index and raw index. the index of a raw is important
            inputmat[i] = inpData['TPM'][idx[0][0]]
    #print(inputmat)
    existing = np.nonzero(inputmat>-1)[0] #finding values in inputmat that are mire than -1 
    #print(existing)
    dists_inp = np.array([mat['dists'][i][0][0] for i in existing])#reduce amount of indexes  
    backgroundgenes_inp = [backgroundgenes['Triv'][i] for i in existing]#reduce amount of indexes
    #print('dists_inp -------------->', dists_inp)
    exp_inp = [inputmat[i] for i in existing]
    exp_inp = np.array(exp_inp)
    #print(max(exp_inp))
    Pvals_low,Pvals_up,Pvals_inter_low,Pvals_inter_up,Pvals_low_adj,Pvals_up_adj, Pvals_inter_low_adj,Pvals_inter_up_adj = CalculatePvalues(dists_inp, exp_inp)
    s = len(exp_inp)
    #return Pvals_low_adj
    print('is printed after the calcpval ------>', Pvals_low[6])
    
    
    Bool_exp = np.zeros(s)
        
    #print(Pvals_low_adj)
    where_are_NaNs = isnan(Pvals_low_adj)
    print('nans -------->', Pvals_low_adj[where_are_NaNs])
    Pvals_low_adj[where_are_NaNs] = 0
    for i in np.nonzero(Pvals_low_adj>0.1)[0]:
        Bool_exp[i]=1
        print(Bool_exp[i])
    output = {backgroundgenes_inp[gene]: Bool_exp[gene] for gene in range(s)}
    
     # writing data into OutputFile as a column of paris gene_name and 1 or zero
    #print(output)
    with open('test.txt', 'w') as f:
        for key,val in output.items():
            f.write('{}, {}\n'.format(key,val))

    return Bool_exp
    #with open(OutputFile) as f:
        #f.write(output)
    #np.savetxt(OutputFile, output)
    #print(np.nonzero(Pvals_low_adj[0]))
