
def BooleanizeFromFile(InputData, OutputFile ):
    """ function for 
        :InputData: str - 
        :OutputFile: str - the name of a file, where the results would be written 
        :return: -
        """
    import scipy.io #need for loading .mat files 
    import numpy as np #need for calculations 
    mat = scipy.io.loadmat('Thresholds_CoTFcRF.mat') # load data from 'Treshholds_CoTFcRF.mat'. Returns a dict, numeric data is under 
                                                     # the key 'dists'. mat['dists'] is a tensor with the shapes : [0:2136][0][0][0:1][0:999]
    
    dt_1 = np.dtype([('Ens', (str, 16)), ('Triv', (str, 10))], align=True) # special datatype for loading data from 'background_genes.txt' 
    with open('background_genes.txt') as inp: 
        backgroundgenes = np.loadtxt(inp, dtype=dt_1)#loading data from 'background_genes.txt'  to a numpy matrix with pairs of values
    print(InputData, OutputFile) 
    dt_2 = np.dtype([('ensembl_gene_id', (str, 20)), ('hgnc_symbol', (str,10)), ('TPM', int)])#special datatype for loading data from InputData 
    with open(InputData) as InpD:
        inpData = np.genfromtxt(InpD, dtype = dt_2, names = True)#loading input data
    
    inputmat = -1*np.ones((2137,1)) #making an array to get rid of genes, that are not in backgroundgenes
    #maybe, it's worth to write len(inpData) instead of 2137
    
    for i in range(2137):#getting rid of genes, that are not in backgroundgenes array
        
        idx = np.nonzero(inpData['ensembl_gene_id']==backgroundgenes['Ens'][i] )#findibg indexis of genes  in inpData, that are not in backgroundgenes

        if len(idx[0])>0: # here and after: np.nonzero returns two arrays: colmn index and raw index. the index of a raw is important
            inputmat[i] = inpData['TPM'][idx[0][0]]

    existing = np.nonzero(inputmat>-1)[0]#finding values in inputmat that are mire than -1 
    #print(existing)
    dists_inp = np.array([mat['dists'][i][0][0] for i in existing])#reduce amount of indexes  
    backgroundgenes_inp = [backgroundgenes['Triv'][i] for i in existing] #reduce amount of indexes
    #print('dists_inp -------------->', dists_inp)
    exp_inp = [inputmat[i] for i in existing] 
    exp_inp = np.array(exp_inp)
    #print(max(exp_inp))
    #calculating Pvals 
    Pvals_low,Pvals_up,Pvals_inter_low,Pvals_inter_up,Pvals_low_adj,Pvals_up_adj, Pvals_inter_low_adj,Pvals_inter_up_adj = CalculatePvalues(dists_inp, exp_inp)
    s = len(exp_inp)
 
    Bool_exp = np.zeros((s,1)) #initializing an array for output 
        
    #print(Pvals_low_adj)
    for i in np.nonzero(Pvals_low_adj>0.1)[0]:
        Bool_exp[i]=1
    output = {backgroundgenes_inp[gene]: Bool_exp[gene] for gene in range(s)}#making a dict for output data 
    
    # writing data into OutputFile as a column of paris gene_name and 1 or zero
    with open(OutputFile, 'w') as f:
        for key,val in output.items():
            f.write('{}, {}\n'.format(key,val[0]))

   

