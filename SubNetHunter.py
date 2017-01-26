# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 10:18:39 2016
This function tries to add each biomass component to the seed and recomputes
paths to the biomass components. This is odne to ensure that all the biomass
components can be reached from the seed set.

@author: aarthi
"""


from __future__ import division
#import networkx as nx
#import multiprocessing
from networkx import read_gpickle, get_node_attributes
#import pdb
#import cPickle
from pickle import dump
import itertools
import time
#import sys
from generate_partitions import generate_partitions
# import cPickle
from operator import mul
#import os
#from check_repetitions import check_repetitions
#from collections import defaultdict
from forward_pass import forward_pass
#import gzip
#import pickle
import math
#from joblib import Parallel, delayed

tic = time.clock()

def calculate_biomass_components(G, src, tar1, alwaysexist1, cutoff, filenames):#, gname2):
    #tic = time.clock()
    #print alwaysexist1
    #print cutoff, 'without 1000 restriction on each'
    H = G.reverse()
    
    nattr = get_node_attributes(G, 'bipartite')
    nattr_inv = {}
    for k, v in nattr.iteritems():
        #print v
        nattr_inv[v] = nattr_inv.get(v, [])
        nattr_inv[v].append(k)  #0 metabolites , 1 are reactions

    for rxns1 in nattr_inv[1]:
        #print rxns1
        if len(set(pred(rxns1)) - alwaysexist1) >= 5:
            G.remove_node(rxns1)
    #pdb.set_trace()
    #if tar1 not in alwaysexist1:
    #print tar1
    tar= tar1
    #pdb.set_trace()
    lowerbound, lowerboundreaction, status_dictfp, me = forward_pass(G, alwaysexist1, src, 0)
    lowerboundrev, lowerboundreactionrev, status_dictrev, me1 = forward_pass(H, me, tar, 1, src)
    status_dict = list(set(status_dictfp).intersection(set(status_dictrev)))
    dagsfound2 = {}
    #pdb.set_trace()
    for item in list(alwaysexist1):
        dagsfound2[item] = {0: ''}
    #cutoff = 25 # sys.argv(1:)

    #print 'cutoff', cutoff

    for rxns in status_dict:
        if set(pred(rxns)).issubset(alwaysexist1):
            for succmets in succ(rxns):
                if succmets not in alwaysexist1:
                    dagsfound2[succmets] = {1: [set([rxns])]}
    #pdb.set_trace()

    for i in range(2, cutoff+1):
        #pdb.set_trace()
        #print i
        
        for rxns in status_dict: #Can you parallalize this?
            #if i == 25 and 'S1 pro_DASH_L_c' in succ(rxns):
             #   pdb.set_trace()
            if set(pred(rxns)).issubset(dagsfound2):

                metsrequired = list(set(pred(rxns)) - alwaysexist1) # This variable avoids the generation of partitions pertaining to seed metabolites
                shortestpathlist = []
                if metsrequired:# and goahead == 'Yes':
                    for predmets in metsrequired:
                        shortestpathlist.append(min(lowerbound[predmets]))
    ##==============================================================================
    ## First loop
    ##==============================================================================
                    for val in range(i-1, len(metsrequired)*(i-2)+1): #range does not include the end value
                        temp = int(math.floor(val/((i-1))))
                        for idx in range(1, temp+1):
                            for combivalue in itertools.combinations(metsrequired, idx):
                                variablemets = list(set(metsrequired) - set(combivalue))
                                temprxnlist = []
                                templen = {}

                                for metabolites in list(combivalue):
                                    if metabolites in dagsfound2:
                                        if i-1 in dagsfound2[metabolites]:
                                            temprxnlist.append([[rxns]])
                                            templen[metabolites]=len(dagsfound2[metabolites][i-1])
                                    else:
                                        break
                                flag = ''
                                if templen:
                                    for mets in list(variablemets):
                                        templen[mets] = len(dagsfound2[mets])
                                    if all([reduce(mul,templen.values()) > 1000, set(succ(rxns)).issubset(set(dagsfound2))]): #, set(succ(rxns)).isdisjoint(set(tar))
                                        flag = 'NA'
                                        #print templen, 'Loop1'
                                    else:
                                        for mets in combivalue:
                                            if i-1 in dagsfound2[mets]:
                                                temprxnlist.append(dagsfound2[mets][i-1])
                                            else:
                                                flag = 'NA'
                                                break
                                    if flag != 'NA':
                                        tempshortestpath = []
                                        for varmet in list(variablemets):
                                            tempshortestpath.append(min(lowerbound[varmet]))

                                        par = generate_partitions(i, tempshortestpath, val-((i-1)*idx))

                                        for partitions in par:
                                            #print partitions,mets
                                            #num += 1
                                            counter = 0
                                            for varmetidx in range(len(variablemets)):
                                                if variablemets[varmetidx] in dagsfound2:
                                                    if partitions[varmetidx] in dagsfound2[variablemets[varmetidx]]:
                                                        counter += 1
                                                        templen[variablemets[varmetidx]] = len(dagsfound2[variablemets[varmetidx]][partitions[varmetidx]])
                                            if counter == len(variablemets): #because we ae appending the value to templen, we need another flag to check if the other metabolite is present or not
                                                if all([reduce(mul,templen.values()) > 1000, set(succ(rxns)).issubset(set(dagsfound2))]): #, set(succ(rxns)).isdisjoint(set(tar))
                                                    break #or pass
                                                else:
                                                    temprxnlist1 = temprxnlist[:]
                                                    for varmetidx in range(len(variablemets)):
                                                        temprxnlist1.append(dagsfound2[variablemets[varmetidx]][partitions[varmetidx]])
                                                    for prodvec in itertools.product(*temprxnlist1):
                                                        rxncomb = set([])
                                                        reactants = set([]) #This is to check for cycles
                                                        for rxnentry in prodvec:
                                                            for individualele in rxnentry:
                                                                if individualele not in rxncomb:
                                                                    rxncomb.add(individualele)
                                                                    for metsreq in pred(individualele):
                                                                        reactants.add(metsreq)
                                                        for succmets in succ(rxns):
                                                            if succmets not in alwaysexist1 and succmets not in reactants:
                                                                if succmets in dagsfound2:
                                                                    if len(rxncomb) in dagsfound2[succmets]:
                                                                        try:
                                                                            if dagsfound2[succmets][len(rxncomb)].index(rxncomb):
                                                                                pass
                                                                        except:
                                                                                dagsfound2[succmets][len(rxncomb)].append(rxncomb)
                                                                    else:
                                                                        dagsfound2[succmets].update({len(rxncomb): [rxncomb]})
                                                                else:
                                                                    dagsfound2[succmets] = {len(rxncomb): [rxncomb]}
    #==============================================================================
    # Second loop
    #==============================================================================
                    for val1 in range(len(metsrequired)*(i-2)+1, len(metsrequired)*(i-1)+1):
                        par1 = generate_partitions(i, shortestpathlist,val1)
                        for partitions in par1:
                                #print partitions,metsrequired
                                temprxnlist2 = []
                                templen = {}
                                flag1 = ''
                                for item in range(len(metsrequired)):
                                    if metsrequired[item] in dagsfound2:
                                        if partitions[item] in dagsfound2[metsrequired[item]]:
                                            templen[metsrequired[item]]=len(dagsfound2[metsrequired[item]][partitions[item]])
                                if templen:
                                    if all([reduce(mul,templen.values()) > 1000, set(succ(rxns)).issubset(set(dagsfound2))]): #, set(succ(rxns)).isdisjoint(set(tar))
                                        flag1 = 'NA'
                                        #print templen
                                    else:
                                        for item in range(len(metsrequired)):
                                            if partitions[item] in dagsfound2[metsrequired[item]]:
                                                temprxnlist2.append(dagsfound2[metsrequired[item]][partitions[item]])
                                            else:
                                                flag1 = 'NA'
                                    if flag1!='NA':
                                        temprxnlist2.append([[rxns]])
                                        for prodvec in itertools.product(*temprxnlist2):
                                            rxncomb1 = set([])
                                            reactants1 = set([]) #This is to check for cycles
                                            for rxnentry in prodvec:
                                                for individualele in rxnentry:
                                                    if individualele not in rxncomb1:
                                                        rxncomb1.add(individualele)
                                                        for metsreq in pred(individualele):
                                                            reactants1.add(metsreq)
                                            for succmets in succ(rxns):
                                                if succmets not in alwaysexist1:# and succmets not in reactants1:
                                                    if succmets in dagsfound2:
                                                        #pdb.set_trace()
                                                        if len(rxncomb1) in dagsfound2[succmets]:
                                                                try:
                                                                    if dagsfound2[succmets][len(rxncomb1)].index(rxncomb1):
                                                                        pass
                                                                except:
                                                                    dagsfound2[succmets][len(rxncomb1)].append(rxncomb1)
                                                        else:
                                                            dagsfound2[succmets].update({len(rxncomb1): [rxncomb1]})
                                                    else:
                                                        dagsfound2[succmets] = {len(rxncomb1): [rxncomb1]}

        #print 'Column value', i, 'Cumulative number of metabolites', len(dagsfound2.keys())
        #toc = time.clock()
        #print 'Time elapsed', toc-tic
        #timetaken = toc-tic
        #counter1 = []
        #if tar1 in dagsfound2:
            #for i in dagsfound2[tar1]:
            #    counter1.append(len(dagsfound2[tar1][i]))
            #print tar, '\t', sum(counter1)
    #for sourcenodes in src:

#        if tar1 in dagsfound2:
##            print tar1, lowerbound[tar1]
#            tarname = os.path.join(filenames[0]+'_'+filenames[1],"".join(src).replace(" ","") +"".join(tar1.split())  +'.pickle.gz')
#            picfname = gzip.open(tarname, 'wb')
#            cPickle.dump([dagsfound2[tar1],cutoff, alwaysexist1],  picfname)
#            picfname.close()
#    with open(gname2[0], 'r') as f:
#        namemap = pickle.load(f)
    #originalrxndags = {}
#    for metabolite in dagsfound2:
#        originalrxndags[metabolite] ={}
#        tempdict1 = {}
#        for plen in dagsfound2[metabolite]:
#            tempdict1[plen] = []
#            #originalrxndags[metabolite][plen] =[]
#            namechanged = []
#            for everyrxnset in dagsfound2[metabolite][plen]:
#                temprxnlist = set([])
#                for element in list(everyrxnset):
#                    temprxnlist.add(namemap[element])
#                namechanged.append(temprxnlist)
#            tempdict1[plen]=namechanged
            
        #originalrxndags[metabolite].update(tempdict1)


    toc = time.clock()
    return lowerbound, lowerboundreaction, dagsfound2, toc-tic
    #,  originalrxndags,
#gname1 = ['S1_P1.gpickle']
#gname1 = sys.argv[1:]
##gname1 = ['/media/aarthi/pras2/GraphsSegregated/BMID000000140386/BMID000000140386_BMID000000140522.gpickle.gz']
cutoff1 = 25
#gname1 = ['/home/aarthi/Dropbox/Consortia_Manuscript/ConsortiaExamples/Simulations/CommunitySImulations/OralCommunity/TwoOrgGraphs/BMID000000140226_BMID000000140256.gpickle']
#gname1 = ['/home/aarthi/Dropbox/Consortia_Manuscript/ConsortiaExamples/Simulations/CommunitySImulations/OralCommunity/SingleGraphs/BMID000000140226.gpickle']
#gname1 = ['/home/aarthi/Dropbox/Consortia_Manuscript/ConsortiaExamples/Simulations/CommunitySImulations/OralCommunity/ShortestPathComparison/140226/BMID000000140226.gpickle']
#gname1 = ['/home/aarthi/Dropbox/Consortia_Manuscript/ConsortiaExamples/Simulations/CommunitySImulations/OralCommunity/ShortestPathComparison/140256/BMID000000140256.gpickle']
#gname1 = ['/home/aarthi/Dropbox/Consortia_Manuscript/ConsortiaExamples/Simulations/CommunitySImulations/PichiaSaccShortestPath/S1_P1_.gpickle']
#gname2 = ['/home/aarthi/Dropbox/Consortia_Manuscript/ConsortiaExamples/Simulations/CommunitySImulations/PichiaSaccShortestPath/S1_P1_namemap.pickle']
#gname1 = ['S1_P1_.gpickle']
#gname2 = ['S1_P1_namemap.pickle']
#gname1 = ['C:\\Users\\Aarthi Ravikrishnan\\Dropbox\\Consortia_Manuscript\\ConsortiaExamples\\Simulations\\PathFindingSimulations\\Cycles\\cycleexample.gpickle']
#gname1 = ['/home/aarthi/Dropbox/Consortia_Manuscript/ConsortiaExamples/Simulations/CommunitySImulations/OralCommunity/ShortestPathComparison/Double/BMID000000140256_BMID000000140226.gpickle']
#gname1 =['C:\Users\Aarthi Ravikrishnan\Dropbox\Consortia_Manuscript\ConsortiaExamples\Simulations\PathFindingSimulations\Examples\PiciBB814\P1.gpickle']
gname1 = ['/home/aarthi/Dropbox/Consortia_Manuscript/ConsortiaExamples/Simulations/PathFindingSimulations/Examples/PiciBB814/P1.gpickle']
G = read_gpickle(gname1[0])
filenames = gname1[0].split('/')[-1].split('.')[0].split('_')
#if not os.path.exists(filenames[0]+'_'+filenames[1]):
#    os.mkdir(filenames[0]+'_'+filenames[1])
#src = [''MNXM99_i','MNXM105_i]
src = ['glc_DASH_D_e']
#src = []
#biomasscomp =  ['S1 gly_c']
#src = ['M1']
#alwaysexist1 = set(['A1','M1'])
#alwaysexist1 = set(['glc_D_e','Ecoli_core_model ACP_c', 'Ecoli_core_model adp_c', 'Ecoli_core_model amp_c', 'Ecoli_core_model atp_c', 'Ecoli_core_model h2o_c', 'Ecoli_core_model h2o_p', 'Ecoli_core_model h_c', 'Ecoli_core_model h_p', 'Ecoli_core_model nad_c', 'Ecoli_core_model nadh_c', 'Ecoli_core_model nadp_c', 'Ecoli_core_model nadph_c', 'Ecoli_core_model pep_c', 'Ecoli_core_model pi_c', 'Ecoli_core_model pi_p', 'Ecoli_core_model ppi_c'])
#alwaysexist1 = set(['Lactococcus_lactis_MG1363 ACP_c', 'Lactococcus_lactis_MG1363 adp_c', 'Lactococcus_lactis_MG1363 amp_c', 'Lactococcus_lactis_MG1363 atp_c', 'Lactococcus_lactis_MG1363 co2_c', 'Lactococcus_lactis_MG1363 coa_c', 'Lactococcus_lactis_MG1363 h2o_c', 'Lactococcus_lactis_MG1363 h2o_p', 'Lactococcus_lactis_MG1363 h_c', 'Lactococcus_lactis_MG1363 h_p', 'Lactococcus_lactis_MG1363 nad_c', 'Lactococcus_lactis_MG1363 nadh_c', 'Lactococcus_lactis_MG1363 nadp_c', 'Lactococcus_lactis_MG1363 nadph_c', 'Lactococcus_lactis_MG1363 pep_c', 'Lactococcus_lactis_MG1363 pi_c', 'Lactococcus_lactis_MG1363 pi_p', 'Lactococcus_lactis_MG1363 ppi_c'])
alwaysexist1 = set(['P1 ACP_c', 'P1 adp_c', 'P1 amp_c', 'P1 atp_c', 'P1 co2_c', 'P1 coa_c', 'P1 h2o_c', 'P1 h2o_p', 'P1 h_c', 'P1 h_p', 'P1 nad_c', 'P1 nadh_c', 'P1 nadp_c', 'P1 nadph_c', 'P1 pep_c', 'P1 pi_c', 'P1 pi_p', 'P1 ppi_c'])
#alwaysexist1 = set(['S1 cmp_c', 'P1 gdp_c', 'P1 zn2_c', 'S1 sel_c', 'P1 fe2_c', 'S1 atp_g', 'S1 dgmp_c', 'P1 h2s_c', 'S1 nh4_c', 'P1 mobd_c', 'S1 na1_c', 'P1 nadp_m', 'P1 atp_c', 'P1 atp_e', 'P1 so4_c', 'S1 nadh_c', 'P1 nadp_c', 'S1 h_v', 'S1 atp_n', 'P1 fad_c', 'P1 fad_m', 'P1 mg2_c', 'S1 h2s_c', 'S1 nadh_r', 'S1 pi_r', 'S1 atp_r', 'P1 gtp_m', 'P1 ump_e', 'S1 atp_x', 'P1 ni2_c', 'S1 pi_c', 'S1 atp_c', 'P1 cu2_c', 'P1 amp_c', 'S1 pi_m', 'P1 coa_r', 'P1 gmp_c', 'S1 dolp_c', 'S1 cbl1_c', 'S1 btn_c', 'P1 k_e', 'P1 k_c', 'P1 coa_c', 'S1 nadp_c', 'P1 nadh_c', 'P1 nadh_m', 'S1 thmpp_c', 'P1 coa_m', 'S1 mobd_c', 'S1 adp_c', 'S1 adp_g', 'P1 4fe4s_c', 'S1 adp_n', 'P1 pi_c', 'P1 adp_c', 'S1 adp_v', 'S1 adp_x', 'P1 adp_m', 'P1 pi_r', 'P1 nadph_r', 'S1 bmocogdp_c', 'S1 10fthf_c', 'P1 cobalt2_c', 'S1 dcmp_c', 'S1 so4_c', 'S1 flxso_c', 'P1 ca2_c', 'S1 atp_m', 'S1 ppi_c', 'P1 nad_m', 'S1 h_m', 'S1 2fe2s_c', 'S1 o2_r', 'S1 nadh_m', 'P1 nad_c', 'P1 co2_m', 'P1 nadph_e', 'P1 nadph_c', 'S1 utp_c', 'S1 o2_c', 'P1 dtmp_c', 'S1 hco3_c', 'P1 co2_c', 'S1 k_e', 'P1 na1_c', 'S1 fe3_c', 'S1 k_c', 'P1 h2o_c', 'P1 ppi_c', 'P1 cl_c', 'S1 nadph_r', 'S1 nad_r', 'P1 pi_m', 'P1 atp_m', 'S1 cl_c', 'S1 nadph_c', 'S1 nad_c', 'P1 btn_c', 'S1 nad_m', 'S1 nadph_m', 'P1 utp_c', 'S1 cobalt2_c', 'P1 nadh_r', 'S1 ca2_c', 'S1 trnaglu_c', 'S1 cu2_c', 'P1 nadph_m', 'P1 2fe2s_c', 'P1 nad_r', 'S1 ctp_c', 'P1 mn2_c', 'P1 thmpp_c', 'S1 adp_m', 'S1 3mop_c', 'S1 4fe4s_c', 'S1 h_x', 'P1 fadh2_m', 'S1 so3_c', 'S1 mg2_c', 'S1 h_r', 'P1 fadh2_c', 'S1 fe2_c', 'P1 fe3_c', 'P1 ctp_c', 'S1 h_n', 'P1 cbl1_c', 'S1 h_c', 'S1 h_e', 'P1 ump_c', 'P1 trnaglu_c', 'P1 damp_c', 'S1 amp_c', 'P1 cdpdag_c', 'S1 co2_m', 'S1 atp_v', 'P1 so3_c', 'P1 nh4_c', 'S1 co2_c', 'S1 mn2_c', 'P1 gtp_c', 'P1 tungs_c', 'S1 zn2_c', 'S1 tungs_c', 'P1 slnt_c', 'S1 coa_r', 'P1 dgmp_c', 'P1 o2_r', 'P1 h_e', 'S1 slnt_c', 'S1 coa_c', 'S1 dtmp_c', 'P1 gdp_m', 'P1 dcmp_c', 'P1 o2_c', 'S1 coa_m', 'P1 h_r', 'P1 sel_c', 'P1 h_v', 'S1 ni2_c', 'P1 h_c', 'P1 bmocogdp_c', 'S1 h2o_c', 'S1 h2o_m', 'P1 dolp_c', 'S1 damp_c', 'P1 h_m'])
#alwaysexistdec = ['bigg_h_i', 'bigg_h2o_i', 'bigg_atp_i', 'bigg_pi_i',
#                  'bigg_nad_i', 'bigg_nadh_i', 'bigg_nadp_i', 'bigg_co2_i',
#                  'bigg_k_i', 'bigg_na1_i','bigg_nadph_i', 'bigg_adp_i',
#                  'bigg_ppi_i', 'bigg_co2_i', 'bigg_coa_i', 'bigg_amp_i',
#                  'bigg_amp_bm', 'bigg_nh3_i', 'bigg_o2_i', 'bigg_glu_L_i',
#                  'bigg_gtp_i', 'bigg_ump_i', 'bigg_udp_i',
#                  'bigg_gdp_i', 'bigg_cmp_i', 'bigg_ctp_i', 'bigg_utp_i',
#                  'bigg_fad_i', 'bigg_gmp_i', 'bigg_dna_i', 'bigg_cdp_i',
#                  'bigg_so3_i', 'bigg_dgtp_i', 'MNXM96063_i', 'bigg_no2_i', 'bigg_h2s_i',
#                  'bigg_imp_i', 'bigg_dadp_i', 'bigg_dgdp_i', 'bigg_itp_i',
#                  'bigg_dcdp_i', 'bigg_datp_i', 'bigg_dctp_i', 'bigg_dtmp_i','bigg_dtmp_bm'
#                  'bigg_dutp_i', 'bigg_dudp_i', 'bigg_cl_i', 'bigg_dgmp_i','bigg_gmp_bm'
#                  'bigg_idp_i', 'bigg_dcmp_i', 'bigg_dump_i', 'bigg_so4_i','bigg_dgmp_bm'
#                  'bigg_pppi_i', 'bigg_damp_i', 'bigg_fmnh2_i', 'MNXM537_i',
#                  'bigg_fe2_i', 'bigg_h2_i', 'bigg_ditp_i', 'bigg_ca2_i',
#                  'bigg_cobalt2_i', 'metacyc_NAD_P_H_i', 'MNXM24_i',
#                  'MNXM90227_i', 'bigg_fadh2_i', 'bigg_nac_i',
#                  'bigg_ncam_i', 'bigg_damp_bm','bigg_dtmp_bm',
#                  'bigg_ump_bm', 'bigg_dcmp_bm', 'bigg_atp_bm', 'metacyc_NAD_P__i',
#                  'bigg_NPmehis_i','MNXM2623_i',
#                  'MNXM11736_i','MNXM11739_i','MNXM96096_i',
#                  'MNXM101_i','bigg_3mop_i','MNXM24_i']
#for alwaysitem in alwaysexistdec:
#    alwaysexist1.add(filenames[0] + ' ' + alwaysitem)
#    alwaysexist1.add(filenames[1] + ' ' + alwaysitem)
G = read_gpickle(gname1[0])
print 'Cutoff', cutoff1
#lowerbound, lowerboundreaction, status_dictfp, me = forward_pass(G, alwaysexist1, src, 0)
#biomass_lb = {}# 'bigg_cmp_bm',
#biomass_components = ['bigg_asn_L_bm', 'bigg_gly_bm', 'bigg_ala_L_bm', 'MNXM18_bm', 'bigg_asp_L_bm', 'bigg_phe_L_bm', 'bigg_ser_L_bm', 'bigg_gln_L_bm', 'bigg_val_L_bm', 'bigg_tyr_L_bm', 'bigg_pro_L_bm',  'MNXM876_bm',  'bigg_trp_L_bm', 'bigg_his_L_bm', 'bigg_ile_L_bm', 'bigg_leu_L_bm',  'bigg_cys_L_bm',  'bigg_thr_L_bm', 'bigg_arg_L_bm', 'bigg_lys_L_bm', 'bigg_met_L_bm' ]
#biomass_components = ['bigg_asn_L_bm', 'bigg_gly_bm', 'bigg_ala_L_bm','MNXM18_bm']
#res_txtfile = os.path.join(filenames[0]+'_'+filenames[1], filenames[0]+'_'+filenames[1]+'_results.txt')
#biomasscomp = []
#for mets in biomass_components:
#    biomasscomp.append(filenames[0]+' ' +mets)
#    biomasscomp.append(filenames[1]+' ' +mets)
#biomasscomp1 = []
#notfound = []
#for item in biomasscomp:
#    #print item
#    if item in lowerbound:
#        biomass_lb[item] = min(lowerbound[item])
#    else:
#        print item, 'not found'
#        notfound.append(item)
#        with open(res_txtfile, 'a') as f:
#            print>>f, item, '\t', 'Not found'
#        biomasscomp.remove(item)
#for mets in notfound:
#    biomasscomp.remove(mets)
#    with open(res_txtfile, 'a') as f:
#        print>>f, item, '\t', 'Not found'

#pdb.set_trace()
#biomass_lb_rev = {}
#for k,v in biomass_lb.iteritems():
#    biomass_lb_rev[v] = biomass_lb.get(v,[])
#for k,v in biomass_lb.iteritems():
#    biomass_lb_rev[v].append(k)
#ikey= min(biomass_lb_rev)
#for biomets in biomass_lb_rev[ikey]:
#    biomasscomp1.append(biomets)
#for startingmets in src:
#    alwaysexist1.add(startingmets)
taritems = ['P1 pyr_c']
#pdb.set_trace()
#res_txtfile = os.path.join(filenames[0]+'_'+filenames[1], filenames[0]+'_'+filenames[1]+'_results.txt')
#for taritems in biomasscomp1:
succ = G.successors
pred = G.predecessors
lb, lbr, dagsfound, timetaken = calculate_biomass_components(G, src, taritems, alwaysexist1, cutoff1, filenames)#, gname2)
print 'Time taken', timetaken
if taritems[0] in dagsfound:
    print len(dagsfound[taritems[0]][cutoff1])
    for i,paths in enumerate(dagsfound[taritems[0]][cutoff1]):
        reactants = set([])
        products = set([])
        for rxns in list(paths):
            for predmets in pred(rxns):
                if predmets not in alwaysexist1:
                   reactants.add(predmets)
            for succmets in succ(rxns):
                if succmets not in alwaysexist1:
                    products.add(succmets)
        if not (reactants-products).issubset(alwaysexist1):# or (products-reactants):#.issubset(alwaysexist1):
                #print i, 'is erroneous'
            print (reactants - products), (products-reactants)
fname = 'pyr' + str(cutoff1)+'.pickle'
with open(fname,'w') as f:
    dump(dagsfound,f)

#while biomasscomp1:
#    for taritems in biomasscomp1:
#        lb, lbr, dagsfound = calculate_biomass_components(G, src, taritems, alwaysexist1, cutoff1, filenames)
#        counter_numberofpaths = []
#        for idx in range(1,cutoff1+1):# dagsfound[taritems]:
#            if idx in dagsfound[taritems]:
#                counter_numberofpaths.append(len(dagsfound[taritems][idx]))
#        with open(res_txtfile, 'a') as f:
#            print>>f, taritems, '\t', sum(counter_numberofpaths)
#        biomasscomp1.remove(taritems)
#        print 'before', len(biomasscomp)
#        biomasscomp.remove(taritems)
#        print 'after', len(biomasscomp)
#        alwaysexist1.add(taritems)
#        if biomasscomp:
#            biomass_lb = {}
#            for item in biomasscomp:
#
#                if item in lb:
#                    biomass_lb[item] = min(lb[item])
#                elif len(biomasscomp) == 1:
#                    biomasscomp.remove(item)
#                else:
#                    print item,'cannot be found'
#            biomass_lb_rev = {}
#            for k,v in biomass_lb.iteritems():
#                biomass_lb_rev[v] = biomass_lb.get(v, [])
#            for k, v in biomass_lb.iteritems():
#                biomass_lb_rev[v].append(k)
#            if biomass_lb_rev:
#                ikey= min(biomass_lb_rev)
#                for biomets in biomass_lb_rev[ikey]:
#                    if biomets not in biomasscomp1:
#                        biomasscomp1.append(biomets)
