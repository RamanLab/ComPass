'''Use the new version of the sbml file for reading the model files from
COBRApy
Modification done - TO remove the reactions where the reactants and the
products are the same from the model itself
COBRApy removes the same reactant and the product , for example
if A+B -> A+C, cobra just gives B->C (which leads to a lot of incorrect
results)
'''
#!/usr/bin/python
#import networkx as networkx
from RevCob import FindExcRxns
import cobra
#import warnings
#warnings.filterwarnings("ignore")
#import collections
#import sys
import os
import glob
import copy
import pdb
#import numpy as np
#import scipy
#from string import ascii_uppercase
organismsdata={}

namemap = {}
#wd = os.getcwd()
def CreateReactions(pathname):#fnames):#,fnames):#fnames can be a list of list of the filenames that have to be considered. If you want to include the files in a folder, use glob
    os.chdir(pathname)
    #pdb.set_trace()
    fnames= glob.glob('*.xml') #pathname ending with /
    organismsdata={}
    print fnames
    #for mn in range(len(fnames)):
    for m in range(len(fnames)):
        [model,useless] = cobra.io.read_sbml_model(fnames[m])#[mn][m]) #To prevent random assignment of dictionary keys and values
        #model, ne = cobra.io.read_sbml_model('/media/aarthi/pras2/Path2MOdels_xml/BMID000000140414.xml')#ecoli_core_model.xml')
        for nonbalrxn in list(set(useless)):
            model.reactions.remove(nonbalrxn)
        temp = cobra.core.ArrayBasedModel(model)
            #org = fnames[m].split('.',1)[0]
        org = model.id
        orgclasdata={org:{'exchange':[],'irrevlhsnodes':[],'irrevrhsnodes':[],'revrhsnodes':[],'revlhsnodes':[],'irrrxnno':[],'revrxnno':[],'exchno':[],'totalnodes':[],'modelrxns':[],'metabolites':[],'exchname':[],'irrevname':[],'revname':[]}}
        #allnamemap = {org:{}}
        #print('1')
        metstmp= []
        rxnids =[]
        mets =[]
        rxnstmp =[]
        #namemap  = {}
            #mets = model.metabolites #This does not work because while passing it, it passes an object
        #print('2')
            #metabch = []
            #mets = model.metabolites
    #==============================================================================
        for metnames in model.metabolites:
            metstmp.append(metnames.id)

        for i in range(len(metstmp)):
            mets.append(metstmp[i]) #In fact, you can pass the name changed metabolites to the function and ask it to work on them instead of changing them later, but this will be a problem with the exchange nodes, so don't do it
    #==============================================================================

        for rxns in model.reactions: #model.reactions is a class of cobra Dictlist
            rxnstmp.append(rxns)

        for i in range(len(rxnstmp)):
            rxnids.append(rxnstmp[i])

        stoi = temp.S
            #lb = temp.lower_bounds
            #ub = temp.upper_bounds
        StoiMat = stoi.toarray()
        StoiMat = stoi.T
        stoi_copy = copy.deepcopy(StoiMat)
        x = StoiMat.toarray()
       # print('3')
            #try:
        revrxn,irrevrxn, excrxn,exc1,irr1,irr2,rev1,rev2,excname, irrevname, revname = FindExcRxns(x,model,temp,mets)
      #  print('4')
        orgclasdata[org]['exchange'] = exc1
        orgclasdata[org]['irrevlhsnodes'] = irr1
        orgclasdata[org]['irrevrhsnodes'] = irr2
        orgclasdata[org]['revlhsnodes'] = rev1
        orgclasdata[org]['revrhsnodes'] = rev2

        orgclasdata[org]['exchname'] = excname
        orgclasdata[org]['irrevname'] = irrevname
        orgclasdata[org]['revname'] = revname
        irrevrxno = []
        for num in range(len(irr1)):
            rno='Org_%s IR' %org + str(num+1)
            irrevrxno.append(rno)
            namemap[rno]= irrevname[num]
        #
        orgclasdata[org]['irrrxnno'] = irrevrxno
        revrxno = []
        for num in range(len(rev1)):
            rno='Org_%s RR' %org + str(num+1)
            revrxno.append(rno)
            namemap[rno]= revname[num]
        orgclasdata[org]['revrxnno'] = revrxno
        revbacrxno = []
        for num in range(len(rev1)):
            rno='Org_%s RevBR' %org + str(num+1)
            revbacrxno.append(rno)
            namemap[rno]= revname[num]

            #excno = []

            #for num in range(len(exc1)):
            #    rno='Org_%s ER' %org + str(num+1)
            #    excno.append(rno)


            #orgclasdata[org]['exchno'] = excno
        #allnamemap[org].update(namemap)
        orgclasdata[org]['totalnodes'] = len(exc1)+len(irr1)+len(rev1)
        orgclasdata[org]['modelrxns']=rxnids
        orgclasdata[org]['revbacrxno'] = revbacrxno
        orgclasdata[org]['metabolites'] = mets
        organismsdata.update(orgclasdata)
        # except:
            #     pass
        #os.chdir(wd)

    return organismsdata, namemap
