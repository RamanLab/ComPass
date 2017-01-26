# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 09:58:43 2016

@author: aarthi
"""
#from joblib import Parallel, delayed
#import networkx as nx
#import math
#import os
import time
from collections import deque, defaultdict
#import cPickle
#import gzip
import pdb

def forward_pass(G, alwaysexist2, src, countertoadd, *args): #tar,
    #print src
   #tic = time.clock()
    targets1 = set([])
    pred = G.predecessors
    succ = G.successors
    alwaysexist = alwaysexist2.copy()


#    alwaysexist = set([])
#    alwaysexist = set(['P1 gdp_c', 'P1 zn2_c', 'S1 sel_c', 'P1 fe2_c', 'S1 atp_g',
#                   'S1 dgmp_c', 'P1 h2s_c', 'S1 nh4_c', 'P1 mobd_c', 'S1 na1_c',
#                   'P1 nadp_m', 'P1 atp_c', 'P1 atp_e', 'P1 so4_c', 'S1 nadh_c',
#                   'P1 nadp_c', 'S1 h_v', 'S1 atp_n', 'P1 fad_c', 'P1 fad_m',
#                   'P1 mg2_c', 'S1 h2s_c', 'S1 nadh_r', 'S1 pi_r', 'S1 atp_r',
#                   'S1 atp_v', 'P1 ump_e', 'S1 atp_x', 'P1 ni2_c', 'S1 pi_c',
#                   'S1 atp_c', 'P1 cu2_c', 'P1 amp_c', 'S1 pi_m', 'P1 coa_r',
#                   'P1 gmp_c', 'S1 dolp_c', 'S1 cbl1_c', 'S1 btn_c', 'P1 k_e',
#                   'P1 k_c', 'P1 coa_c', 'S1 nadp_c', 'P1 nadh_c', 'P1 nadh_m',
#                   'S1 thmpp_c', 'P1 coa_m', 'S1 mobd_c', 'S1 adp_c', 'P1 4fe4s_c',
#                   'P1 pi_c', 'P1 adp_c', 'P1 pi_r', 'P1 nadph_r', 'S1 bmocogdp_c',
#                   'S1 10fthf_c', 'P1 cobalt2_c', 'S1 dcmp_c', 'S1 so4_c',
#                   'S1 flxso_c', 'P1 ca2_c', 'S1 atp_m', 'S1 ppi_c', 'P1 nad_m',
#                   'S1 h_m', 'S1 2fe2s_c', 'S1 o2_r', 'S1 nadh_m', 'P1 nad_c',
#                   'P1 co2_m', 'P1 nadph_e', 'P1 nadph_c', 'S1 utp_c', 'S1 o2_c',
#                   'P1 dtmp_c', 'S1 hco3_c', 'P1 co2_c', 'S1 k_e', 'P1 na1_c',
#                   'S1 fe3_c', 'S1 k_c', 'P1 h2o_c', 'P1 ppi_c', 'P1 cl_c',
#                   'S1 nadph_r', 'S1 nad_r', 'P1 pi_m', 'P1 atp_m', 'S1 cl_c',
#                   'S1 nadph_c', 'S1 nad_c', 'P1 btn_c', 'S1 nad_m', 'S1 nadph_m',
#                   'P1 utp_c', 'S1 cobalt2_c', 'P1 nadh_r', 'S1 ca2_c', 'S1 trnaglu_c',
#                   'S1 cu2_c', 'P1 nadph_m', 'P1 2fe2s_c', 'P1 nad_r', 'S1 ctp_c',
#                   'P1 mn2_c', 'P1 thmpp_c', 'S1 3mop_c', 'S1 4fe4s_c', 'S1 h_x',
#                   'P1 fadh2_m', 'S1 so3_c', 'S1 mg2_c', 'S1 h_r', 'P1 fadh2_c',
#                   'S1 fe2_c', 'P1 fe3_c', 'P1 ctp_c', 'S1 h_n', 'P1 cbl1_c',
#                   'S1 h_c', 'S1 h_e', 'P1 ump_c', 'P1 trnaglu_c', 'P1 damp_c',
#                   'S1 amp_c', 'P1 cdpdag_c', 'S1 co2_m', 'P1 so3_c', 'P1 nh4_c',
#                   'S1 co2_c', 'S1 mn2_c', 'P1 gtp_c', 'P1 tungs_c', 'S1 zn2_c',
#                   'S1 tungs_c', 'P1 slnt_c', 'S1 coa_r', 'P1 dgmp_c', 'P1 o2_r',
#                   'P1 accoa_c', 'P1 h_e', 'S1 slnt_c', 'S1 coa_c', 'S1 dtmp_c',
#                   'P1 dcmp_c', 'P1 o2_c', 'S1 coa_m', 'P1 h_r', 'P1 sel_c',
#                   'P1 h_v', 'S1 ni2_c', 'P1 h_c', 'P1 bmocogdp_c', 'S1 h2o_c',
#                   'S1 h2o_m', 'P1 dolp_c', 'S1 damp_c', 'P1 h_m', 'S1 gtp_c',
#                   'P1 ppi_m', 'P1 h2o_m','P1 s_c','S1 s_c','S1 imp_c','S1 co2_n','S1 h2o_x','S1 adp_m'
#                   ,'S1 coa_x', 'S1 ump_c','S1 amp_x', 'S1 ppi_x','S1 cmp_c','P1 adp_m','S1 h2o_v']) #added 'P1 s_c','S1 imp_c','S1 co2_n'
#
    for ele in args:
        for targets in ele:

            targets1.add(targets)

            if targets in alwaysexist:
                alwaysexist.remove(targets)
    if countertoadd == 1:
        for srcs in src:
            for tarrxns in succ(srcs):
                for tarreqmets in pred(tarrxns):
                    alwaysexist.add(tarreqmets)
                    #print tarreqmets
        #alwaysexist.add
    #print alwaysexist
    lowerbound = defaultdict(list)
    lowerboundreaction = defaultdict(list)
    for item in alwaysexist:
        lowerbound[item].append(0)
    stage = 1
    #if targets not in alwaysexist:
        #visiteddisc, reactions, pathsfinaldict, stuck_dict= (defaultdict(list) for i in range(4))
    queue = deque([])
    status_dict = defaultdict(str)
    #mediacomp = (set([]) for i in range(1))
    alwaysexist.add(src[0])
    mediacomp = alwaysexist.copy()
    startnode = []
    cannotbetriggered = defaultdict(str)
    for startingnodes in alwaysexist:
        if startingnodes in G:
            for startingrxns in succ(startingnodes):

                if set(pred(startingrxns)).issubset(alwaysexist):
                    startnode.append(startingrxns)
                    for metsprod in succ(startingrxns):
                        mediacomp.add(metsprod)
                        if stage not in lowerbound[metsprod]:
                            lowerbound[metsprod].append(stage)
                    if stage not in lowerboundreaction[startingrxns]:
                        lowerboundreaction[startingrxns].append(stage)
    for startrxn in startnode:
        temp = succ(startrxn)
        for item in temp:
            for rxns in succ(item):
                if set(pred(rxns)).issubset(mediacomp):# and set(pred(rxns)).isdisjoint(targets1):
                    queue.append(rxns)
                else:
                    cannotbetriggered[rxns] = 'Y'
        status_dict[startrxn] = 'V'
    while queue:
        stage += 1
        for parent in list(queue):
            #if parent == 'R4':
             #   pdb.set_trace()
            if status_dict[parent] == '':


                temp = succ(parent)
                if stage not in lowerboundreaction[parent]:
                    lowerboundreaction[parent].append(stage)
                for item in temp:
                    mediacomp.add(item)
#                    if item == 'M33':
#                        pdb.set_trace()
                    if stage not in lowerbound[item]:
                        lowerbound[item].append(stage)
                    for progeny in succ(item):
#                        if progeny == 'R16':
#                            pdb.set_trace()
                        if set(pred(progeny)).issubset(mediacomp):# and set(pred(progeny)).isdisjoint(targets1):
                            if status_dict[progeny] != 'V':
                                queue.append(progeny)
                                status_dict[parent] = 'V'
                            #else: #Reomve it to get one cycles
                            #    status_dict[parent] = 'V'

                        else:
                            cannotbetriggered[progeny] = 'Y'
            elif status_dict[parent] == 'V':
                for item2 in succ(parent):
                    if stage not in lowerbound[item2]:
                        lowerbound[item2].append(stage)
            queue.popleft()

    #toc = time.clock()
    return lowerbound, lowerboundreaction, status_dict, mediacomp
    #return mediacomp, toc - tic, alwaysexist, status_dict, lowerbound, lowerboundreaction
#gname1 = ['/home/aarthi/Dropbox/Graphs/PichiaSaccharomyces/S1_P1.gpickle']
#G = nx.read_gpickle(gname1[0])
#me, time1, ae, status_dict, lowerbound, lowerboundreaction = biomass_comp_pathfind(G, ['S1 glc_DASH_D_c'], 'xyltoetoh', 'S1 etoh_c')
#founddic = defaultdict(list)
#with open('/home/aarthi/Desktop/Singleorgnames.txt','r') as f:
#    fnames1 = f.read().splitlines()
##fnames2 = ['BMID000000140586.gpickle']
#for names in fnames1:
#    print names
#    pathname = os.path.join('/media/aarthi/pras2/Path2ModelsSingleGraph1/',names)
#    G = nx.read_gpickle(pathname)
#    filenames = pathname.split('/')[-1].split('.')[0].split('_')
#    biomass_components = ['bigg_asn_L_bm', 'bigg_gly_bm', 'bigg_ala_L_bm', 'MNXM18_bm', 'bigg_asp_L_bm', 'bigg_phe_L_bm', 'bigg_ser_L_bm', 'bigg_gln_L_bm', 'bigg_val_L_bm', 'bigg_tyr_L_bm', 'bigg_pro_L_bm',  'MNXM876_bm',  'bigg_trp_L_bm', 'bigg_his_L_bm', 'bigg_ile_L_bm', 'bigg_leu_L_bm',  'bigg_cys_L_bm',  'bigg_thr_L_bm', 'bigg_arg_L_bm', 'bigg_lys_L_bm', 'bigg_met_L_bm' ]
#    pred = G.predecessors
#    succ = G.successors
#    src = [filenames[0]+ ' ' + 'MNXM99_i']
#    tarmets = []
#    found = []
#    for item in biomass_components:
#        tarmets.append(filenames[0]+ ' ' +item)
#    for tar in tarmets:
#        me, time1, ae = biomass_comp_pathfind(G, src, filenames, tar)
#        if tar in me:
#            found.append(tar)
#    founddic[names] = found
#print time1
#
##with open('/home/aarthi/Dropbox/Graphs/Results/NewimplementationScope/NewimplementationScope1.gpickle','w') as f:
##    cPickle.dump(founddic,f)
#
#import pdb
#from collections import deque, defaultdict
#def forward_pass(G, alwaysexist, src):
#
#    for startingmetnodes in src:
#        alwaysexist.add(startingmetnodes)
#    alwaysexist3 = alwaysexist.copy()
#    succ = G.successors
#    pred = G.predecessors
#    queue = deque([])
#    status_dict = defaultdict(str)
#    cannotbetriggered = defaultdict(str)
#    stage = 1
#    lowerbound = defaultdict(list)
#    lowerboundreaction = defaultdict(list)
#    for startingnodes in list(alwaysexist):
#        if startingnodes == 'Org_BMID000000140621 IR680':
#            pdb.set_trace()
#        if startingnodes in G: #alwaysexist set is a univeral one in biomodels analysis -- so this check is essential
#            for startingrxns in succ(startingnodes):
#                #status_dict[startingrxns] = 'V'
#                if set(pred(startingrxns)).issubset(alwaysexist):
#                    if stage not in lowerboundreaction[startingrxns]:
#                        lowerboundreaction[startingrxns].append(stage)
#                    queue.append(startingrxns)
#                    for metsprod in succ(startingrxns):
#                        alwaysexist3.add(metsprod)
#                        if stage not in lowerbound[metsprod]:
#                            lowerbound[metsprod].append(stage)
#    #pdb.set_trace()
#    while queue:
#        stage += 1
#        #print stage, len(queue)
#        for parent in list(queue):
#            if parent == 'Org_S1 RR377':
#                pdb.set_trace()
#            if status_dict[parent] == '':
#                status_dict[parent] = 'V'
#                for item in succ(parent):
#                    if stage not in lowerbound[item]:
#                        lowerbound[item].append(stage)
#                    for progeny in succ(item):
#                        if set(pred(progeny)).issubset(alwaysexist3):
#                            if status_dict[progeny] != 'V':
#                                queue.append(progeny)
#                            elif status_dict[parent] == 'V':
#
#                                for item in succ(parent):
#                                    if stage not in lowerbound[item]:
#                                        lowerbound[item].append(stage)
#                        else:
#                            cannotbetriggered[progeny] = 'Y'
#                    alwaysexist3.add(item)
#            if stage not in lowerboundreaction[parent]:
#                lowerboundreaction[parent].append(stage)
#            queue.popleft()
#    pdb.set_trace()
#    return lowerbound, lowerboundreaction, status_dict
