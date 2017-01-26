#!/usr/bin/python
import pdb
import networkx as nx
#from networkx.algorithms import bipartite
#import matplotlib.pyplot as plt
import CreateGraphFromS
import itertools
#import itertools
import os
import sys
from cPickle import dump
#n is the number of organisms that need to be in consortia
#cluster them based on reactions and choose the organisms from different clusters
#class GraphCreate(orgs):
excrxnno=[]
excrxnnor=[]
def InternalCreate(organismsdata):
    #combos = list(itertools.combinations(orgs, n))
    orgkeys1 = organismsdata.keys()
    G = nx.DiGraph()
   # pdb.set_trace()
    for j in range(len(orgkeys1)):
    #    rxnkeys = organismsdata[orgkeys[j]].keys()
        G.add_nodes_from(organismsdata[orgkeys1[j]]['exchno'],bipartite=1)
        G.add_nodes_from(organismsdata[orgkeys1[j]]['irrrxnno'],bipartite=1)
        G.add_nodes_from(organismsdata[orgkeys1[j]]['revrxnno'],bipartite = 1)
        G.add_nodes_from(organismsdata[orgkeys1[j]]['revbacrxno'],bipartite = 1)
        irrevlhs = list(set([item for sublist in organismsdata[orgkeys1[j]]['irrevlhsnodes'] for item in sublist]))
        irrevrhs = list(set([item for sublist in organismsdata[orgkeys1[j]]['irrevrhsnodes'] for item in sublist]))
        revlhs = list(set([item for sublist in organismsdata[orgkeys1[j]]['revlhsnodes'] for item in sublist]))
        revrhs = list(set([item for sublist in organismsdata[orgkeys1[j]]['revrhsnodes'] for item in sublist]))
        G.add_nodes_from(irrevlhs,bipartite = 0)
        G.add_nodes_from(irrevrhs,bipartite = 0)
        G.add_nodes_from(revlhs,bipartite = 0)
        G.add_nodes_from(revrhs,bipartite = 0)

        #G.add_nodes_from(organismsdata[orgkeys[j]]['revlhsnodes'],bipartite = 0)
        #G.add_nodes_from(organismsdata[orgkeys[j]]['irrevrhsnodes'],bipartite = 0)
        #G.add_nodes_from(organismsdata[orgkeys[j]]['revrhsnodes'],bipartite = 0)
        ###Irreversible module####

        for idx in range(len(organismsdata[orgkeys1[j]]['irrrxnno'])):
            for idx2 in range(len(organismsdata[orgkeys1[j]]['irrevlhsnodes'][idx])):
                G.add_edges_from([(organismsdata[orgkeys1[j]]['irrevlhsnodes'][idx][idx2],organismsdata[orgkeys1[j]]['irrrxnno'][idx])])
            for idx3 in range(len(organismsdata[orgkeys1[j]]['irrevrhsnodes'][idx])):
                G.add_edges_from([(organismsdata[orgkeys1[j]]['irrrxnno'][idx],organismsdata[orgkeys1[j]]['irrevrhsnodes'][idx][idx3])])

        ###Reversible module####
        for idx in range(len(organismsdata[orgkeys1[j]]['revrxnno'])):
             for idx2 in range(len(organismsdata[orgkeys1[j]]['revlhsnodes'][idx])):
                 G.add_edges_from([(organismsdata[orgkeys1[j]]['revlhsnodes'][idx][idx2],organismsdata[orgkeys1[j]]['revrxnno'][idx])])
                 G.add_edges_from([(organismsdata[orgkeys1[j]]['revbacrxno'][idx],organismsdata[orgkeys1[j]]['revlhsnodes'][idx][idx2])])
             for idx3 in range(len(organismsdata[orgkeys1[j]]['revrhsnodes'][idx])):
                 G.add_edges_from([(organismsdata[orgkeys1[j]]['revrxnno'][idx],organismsdata[orgkeys1[j]]['revrhsnodes'][idx][idx3])])
                 G.add_edges_from([(organismsdata[orgkeys1[j]]['revrhsnodes'][idx][idx3],organismsdata[orgkeys1[j]]['revbacrxno'][idx])])
    return G

def ExchangeCreate(G,orgs,namemap):
#    metstrip = []
    #pdb.set_trace()
    orgkeys1 = orgs.keys()
    mexcids = []
 #   pdb.set_trace()
    for j in range(len(orgkeys1)):
        #me = organismsdata[orgkeys[j]]['metabolites']
        ec = orgs[orgkeys1[j]]['exchange']
        mexcids.append(ec)
    #pdb.set_trace()
    commonexc = list(set.intersection(*map(set,mexcids))) #Common exchange metabolites in different organisms

    for j in range(len(orgkeys1)):
        renamedexc = [orgkeys1[j] + ' ' + s for s in commonexc]
        #pdb.set_trace()
        exclen = range(0,len(commonexc))
        excrxnno = ['Org_%s ER' %orgkeys1[j] +  str(t+1) for t in exclen]
        #excrxntoreturn.append(excrxnno)
        excrxnnor = ['Org_%s ERR' %orgkeys1[j] + str(t+1) for t in exclen]
        #excrxnnortoreturn.append(excrxnnor)
        G.add_nodes_from(excrxnno, bipartite=1)
        G.add_nodes_from(excrxnnor, bipartite=1)
        G.add_nodes_from(commonexc, bipartite=0)
        G.add_nodes_from(renamedexc, bipartite=0)
#        pdb.set_trace()
        for k in range(len(renamedexc)):
            namemap[excrxnno[k]] = commonexc[k]
            namemap[excrxnnor[k]] = commonexc[k]
            G.add_edges_from([(renamedexc[k],excrxnno[k])])
            G.add_edges_from([(excrxnno[k],commonexc[k])])
            G.add_edges_from([(commonexc[k],excrxnnor[k])])
            G.add_edges_from([(excrxnnor[k],renamedexc[k])])
    #pdb.set_trace() #3 orgs and 2 orgs
    for j in range(len(orgkeys1)):
        metitems = orgs[orgkeys1[j]]['exchange']
        noncommonexc = list(set(metitems) - set(commonexc))
        ncrenamedexc = [orgkeys1[j] + ' ' + s for s in noncommonexc]
        ncexclen = range(0,len(noncommonexc))
        ncexcrxnno = ['Org_%s NCER' %orgkeys1[j] +  str(t+1) for t in ncexclen]
        ncexcrxnnor = ['Org_%s NCERR' %orgkeys1[j] + str(t+1) for t in ncexclen]

        G.add_nodes_from(ncexcrxnno, bipartite=1)
        G.add_nodes_from(ncexcrxnnor, bipartite=1)
        G.add_nodes_from(noncommonexc, bipartite=0)
        G.add_nodes_from(ncrenamedexc, bipartite=0)
        for k in range(len(ncrenamedexc)):
            namemap[ncexcrxnno[k]] = noncommonexc[k]
            namemap[ncexcrxnnor[k]] = noncommonexc[k]


            G.add_edges_from([(ncrenamedexc[k],ncexcrxnno[k])])
            G.add_edges_from([(ncexcrxnno[k],noncommonexc[k])])
            G.add_edges_from([(noncommonexc[k],ncexcrxnnor[k])])
            G.add_edges_from([(ncexcrxnnor[k],ncrenamedexc[k])])


    return G, excrxnno, excrxnnor,noncommonexc, ncexclen,ncrenamedexc,commonexc,namemap
#
#
#
#if __name__ == '__main__':
print 'Enter full path'
pname= raw_input()
#pname = '/media/aarthi/pras2/Path2ModelsForgraph/'
organismsdata, allnamemap = CreateGraphFromS.CreateReactions(pname)#,fname)
#print 'Enter the number of organisms to be in consortia'
#kno = raw_input()
kno = 2
orgkeys = organismsdata.keys()

#if kno == 1:
#    combolist = orgkeys
#else:
combolist = list(itertools.combinations(range(len(orgkeys)),int(kno)))
#pdb.set_trace()
if combolist:
    for m in range(len(combolist)):
        fnm = ''
        tempdict= {}
        for n in range(len(combolist[m])):
            tempdict[orgkeys[combolist[m][n]]] = organismsdata[orgkeys[combolist[m][n]]]
            fnm = fnm +  orgkeys[combolist[m][n]] + '_'
        G = InternalCreate(tempdict)
        H, excnos1, excnos2,A,B,C,commonexc,namemap = ExchangeCreate(G,tempdict,allnamemap)
        #fnm ='Combination' + str(m)
        #fnm = fnm + 'windows'
    #        pname1 = '/media/aarthi/pras2/Path2ModelsSingleGraph/'
    #        filename = os.path.join(pname1,fnm)
        nx.write_gpickle(H,fnm+ '.gpickle')
        print 'len edges', len(G.edges())
        print 'len nodes', len(G.nodes())
        with open(fnm+'namemap'+'.pickle', 'w') as f:
            dump(namemap, f)
    sys.path.append(pname)
else:
    print 'Number of organisms for consortia is more than the models given'
#if __name__ == "__main__":
#    InternalCreate(orgs)
#    ExchangeCreate(G,orgs)


        # for j in range(len(orgkeys)):
        #     for item in range(len(organismsdata[orgkeys[j]]['irrevlhsnodes'])):
        #         if organismsdata[orgkeys[j]]['irrevlhsnodes'][item]!=[] and organismsdata[orgkeys[j]]['irrevrhsnodes'][item]!=[]:
        #             #rxnno='Org_%s IR'%orgkeys[j]+str(item+1)

        #             for item2 in range(len(organismsdata[orgkeys[j]]['irrevlhsnodes'][item])):
        #                 G.add_nodes_from(organismsdata[orgkeys[j]]['irrevlhsnodes'][item][item2],bipartite=0)
        #                 G.add_nodes_from(organismsdata[orgkeys[j]]['irrrxnnodes'],bipartite=1)
        #                 G.add_edges_from([(organismsdata[orgkeys[j]]['irrevlhsnodes'][item][item2],organismsdata[orgkeys[j]]['irrrxnnodes'])])
        #             for item3 in range(len(organismsdata[orgkeys[j]]['irrevrhsnodes'][item])):
        #                 G.add_nodes_from([organismsdata[orgkeys[j]]['irrevrhsnodes'][item][item3]],bipartite=0)
        #                 G.add_edges_from([(rxnno, organismsdata[orgkeys[j]]['irrevrhsnodes'][item][item3])])

        #     #Adding reversible nodes
        #     for item4 in range(len(organismsdata[orgkeys[j]]['revrhsnodes'])):
        #         revrxn='Org_%s RR'%orgkeys[j]+str(item4+1)
        #         revrxn1='Org_%s ReR'%orgkeys[j]+str(item4+1)
        #         for item5 in range(len(organismsdata[orgkeys[j]]['revlhsnodes'][item4])):
        #             G.add_nodes_from(organismsdata[orgkeys[j]]['revlhsnodes'][item4][item5],bipartite=0)
        #             G.add_nodes_from(revrxn,bipartite=1)
        #             G.add_nodes_from(revrxn1,bipartite=1)
        #             G.add_edges_from([(organismsdata[orgkeys[j]]['revlhsnodes'][item4][item5],revrxn)])
        #             G.add_edges_from([(revrxn1,organismsdata[orgkeys[j]]['revlhsnodes'][item4][item5])])
        #         for item6 in range(len(organismsdata[orgkeys[j]]['revrhsnodes'][item4])):
        #             G.add_nodes_from(organismsdata[orgkeys[j]]['revrhsnodes'][item4][item6],bipartite=0)
        #             G.add_edges_from([(revrxn,organismsdata[orgkeys[j]]['revrhsnodes'][item4][item6])])
        #             G.add_edges_from([(organismsdata[orgkeys[j]]['revrhsnodes'][item4][item6],revrxn1)])
        # return G
