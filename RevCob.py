#!/usr/bin/python
import cobra
import numpy as np
import scipy
import copy
###FINAL VERSION###
#model, ne = cobra.io.read_sbml_model('/media/aarthi/pras2/Path2MOdels_xml/BMID000000140414.xml')#ecoli_core_model.xml')
#temp = cobra.core.ArrayBasedModel(model)
#for nonbalrxn in list(set(ne)):
#    model.reactions.remove(nonbalrxn)
#metstmp= []
#rxnids =[]
#mets =[]
#rxnstmp =[]
#for i in model.metabolites:
#    metstmp.append(i)
#
#for i in range(len(metstmp)):
#    mets.append(str(metstmp[i]))
#for i in model.reactions:
#    rxnstmp.append(i)
#
#for i in range(len(rxnstmp)):
#     rxnids.append(str(rxnstmp[i]))
#
#stoi = temp.S
#lb = temp.lower_bounds
#ub = temp.upper_bounds
#StoiMat = stoi.toarray()
#StoiMat = stoi.T
#stoi_copy = copy.deepcopy(StoiMat)
#x = StoiMat.toarray()
#Since array operations in python are much simple when done as rows,
#we take the transpose of the stoihiometric matrix
#The S matrix will now be Reactions * Metabolites
#
'''Here x is the stoichiometric matrix
model is the cobra model
temp is the cobra array based model
mets is the model metabolites'''
def FindExcRxns(x,model,temp,mets):
#    print(model.id)
    #x[np.where(x==-1)]==1

    xdim,ydim = np.shape(x)
    Amat=[]
    Bmat=[]
    Cmat=[]
    Amat_len=[]
    Bmat_len=[]
    Cmat_len=[]
    ExcIdx=[]
    rxns=[]
    rxnNames=[]
    for j in model.reactions:
        rxns.append(j.reaction)
    for j in range(len(model.reactions)):
        dumb = model.reactions[j]
        rxnNames.append(dumb.id)

    for i in range(xdim):
        Amat.append(np.where(x[i]==-1))
        Bmat.append(np.where(x[i]!=0))
        Cmat.append(np.where(x[i]==1))
        Amat_len.append(len(Amat[i][0]))
        Bmat_len.append(len(Bmat[i][0]))
        Cmat_len.append(len(Cmat[i][0]))
        #Case 1 - Presence of bulk metabolites in the medium
        if rxns[i][-1] == 'b': #Assuming the bulk metabolites end in 'b'
            if Amat_len[i]==1 and Cmat_len[i]==1:
                ExcIdx.append(i)

        #Case 2 - Presence of exchange metabolites
        elif Amat_len[i]==1 and Bmat_len[i]==1:
            ExcIdx.append(i)
        elif Cmat_len[i]==1 and Bmat_len[i]==1:
            ExcIdx.append(i)
    tempmets = []
    excnodes = []
    excmetids =[]
    exchids=[]
    excnodestemp =[]
    for l in ExcIdx:
        exchids.append(rxnNames[l])
        if rxns[l][-1] == 'b':
            excnodes.append(mets[np.nonzero(x[l])[0][0]]) # This will generate a list of tuples
            #for m in excnodestemp:
            #    excnodes.append(mets[m])
        else:
            excmetids.append(np.nonzero(x[l]))
            #continue
    # for l in ExcIdx:
    #     excmetids.append(np.nonzero(x[l]))
    if excmetids:
        for k in range(len(excmetids)):
            dummy2 =excmetids[k][0].tolist() # To convert it to an array
            excnodestemp.append(dummy2)
    #dummy=[]
    #excnodestemp1 = sum(excnodestemp,[])
        excnodestemp1=[item4 for sublist in excnodestemp for item4 in sublist]
        for m in excnodestemp1:
            excnodes.append(mets[m])
    #tempmets.append(mets[l])
    #excnodes.append(tempmets)
    #tempmets =[]

    #return ExcIdx,rxns

    #def InternalRxns(x):
    AllRxnsIds = []
    for i in range(len(rxns)):
        AllRxnsIds.append(i)

    InternalRxns=list(set(AllRxnsIds)^set(ExcIdx))

    ReversibleRxns = []
    IrreversibleRxns = []
    lb = temp.lower_bounds
    ub = temp.upper_bounds

    for i in InternalRxns: #Changed to InternalRxns from range(len(InternalRxns))
        if lb[i] < 0 and ub[i] >= 0:
            ReversibleRxns.append(i)
        elif lb[i] >=0 and ub[i] >=0:
            IrreversibleRxns.append(i)
    #Irreversible part
    irrevlhstemp1 = []
    irrevlhstemp =[]
    irrevlhsnodes=[]
    tempmets =[]
    irrevrxnids = []
    #irrevlhsnodes1 = []
    for j in IrreversibleRxns:
        irrevrxnids.append(rxnNames[j])
        irrevlhstemp1.append(np.where(x[j]<0))
        #if 335 <= j <= 337:
        #    print(irrevlhstemp1[-5:])
    for k in range(len(irrevlhstemp1)):
        dummy =irrevlhstemp1[k][0].tolist()
        irrevlhstemp.append(dummy)
    for l in range(len(irrevlhstemp)):
        for m in irrevlhstemp[l]:
            metch=model.id + ' ' + mets[m]
            tempmets.append(metch)
        irrevlhsnodes.append(tempmets)
        tempmets =[]
    # for l in range(len(irrevlhsnodes1)):
    #     for m in irrevlhsnodes1[l]:
    #     #for num in irrevlhsnodes1:
    #         metch=model.id + ' ' + irrevlhsnodes1[l][m] #Can the same variable be used again? Can 'i' of the first iteration used in the loop of second? i.e the variable name use can be repeated??
    #         irrevlhsnodes.append(metch)

    tempmets=[]
    irrevrhsnodes=[]
    #irrevrhsnodes1=[]
    irrevrhstemp1 = []
    irrevrhstemp =[]
    for j in IrreversibleRxns:
        irrevrhstemp1.append(np.where(x[j]>0)) #Empty are the ones where the stoichiometric matrix itself has nothing
    for k in range(len(irrevrhstemp1)):
        dummy =irrevrhstemp1[k][0].tolist()
        irrevrhstemp.append(dummy)
    for l in range(len(irrevrhstemp)):
        for m in irrevrhstemp[l]:
            metch=model.id + ' ' + mets[m]
            tempmets.append(metch)
        irrevrhsnodes.append(tempmets)
        tempmets =[]
    # for num in irrevrhsnodes1:
    #     metch=model.id + ' ' + num
    #     irrevrhsnodes.append(metch)

    #Reversible part
    revlhsnodes =[]
    #revlhsnodes1 = []
    revlhstemp1 = []
    revlhstemp =[]
    revrxnids = []
    for j in ReversibleRxns:
        revrxnids.append(rxnNames[j])
        revlhstemp1.append(np.where(x[j]<0))
    for k in range(len(revlhstemp1)):
        dummy =revlhstemp1[k][0].tolist()
        revlhstemp.append(dummy)
        #print('*')
    for l in range(len(revlhstemp)):
        for m in revlhstemp[l]:
            metch=model.id + ' ' + mets[m]
            tempmets.append(metch)
        revlhsnodes.append(tempmets)
        tempmets =[]
    # for num in revlhsnodes1:
    #     metch=model.id + ' ' + num
    #     revlhsnodes.append(metch)
    revrhsnodes =[]
    revrhstemp1 = []
    revrhstemp =[]
    #revrhsnodes1 = []
    for j in ReversibleRxns:
        revrhstemp1.append(np.where(x[j]>0))
    for k in range(len(revrhstemp1)):
        dummy =revrhstemp1[k][0].tolist()
        revrhstemp.append(dummy)
    for l in range(len(revrhstemp)):
        for m in revrhstemp[l]:
            metch=model.id + ' ' + mets[m]
            tempmets.append(metch)
        revrhsnodes.append(tempmets)
        tempmets =[]
        # for num in revrhsnodes1:
        #     metch=model.id + ' ' + num
    return ReversibleRxns, IrreversibleRxns, ExcIdx, excnodes, irrevlhsnodes, irrevrhsnodes, revlhsnodes, revrhsnodes, exchids, irrevrxnids, revrxnids
# for k in range(len(irrevlhstemp[0][0])):
#     for l in irrevlhstemp[0][0][k]:
#         templist.append(mets[l])
#         templist =[]
#     Irrevnodes.append(templist)

#RevIntRxns=list(set(ReversibleRxns)^set(ExcIdx))

    #Amat.append(StoiMat[i][np.nonzero(StoiMat[i]<0)])
    #Amat_sum.append(Amat[i].sum())
    #idx = np.where(StoiMat[i]==-1)
#        if idx:
#            StoiMat[i][idx] = 1
#            AddedStoi1 = len(idx[0])
#            Amat.append(AddedStoi1)
    #if Amat != 0:
    #    for x in range(len(Amat)):
    #        Amat[i]=1Bmat
# Bmat = []
# for j in range(xdim):
#     idx2 = np.where(StoiMat[j]!=0)
#     if idx2[0]:
#         StoiMat[i][idx2] = 1
#         NonZeroStoi1 = len(idx2[0])
#         Bmat.append(NonZeroStoi1)
