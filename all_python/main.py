import qgsm_collect as collect
import qgsm_FB_corr as fb_corr
import numpy as np

energy = 7000
#O = collect.QGSM_Distributions(energy)
#O.collectData(bcorr=True,nbnf=False)
#O.store_bin_data()
#O.create_nbnf_hists()
#O.writeAnalysis()

bcorr = np.loadtxt("bcorr_temp.out",delimiter=',')
####900 GeV
##NeventsList = [785304, 634299, 599483]
##linecount = 2948250
####7000 GeV
##NeventsList = [3016648, 2456139, 2369212]
##linecount = 18251334
FBs = fb_corr.FB(energy)
##FBs.b_corr_count(NeventsList[2],bcorr,linecount)
#FBs.b_corr_collect_file()
#FBs.b_corr_collect_mem(O.bcorr_Nevents,O.bcorr,O.bcorr_linecount)
#FBs.b_corr_count(O.bcorr_Nevents,O.bcorr,O.bcorr_linecount)
FBs.b_corr_plot(from_file='bcorr.py.out')

#nB_nFs = fb_corr.FB(energy)
#nB_nFs.nB_nF_collect_mem(O.nbnf_all_linecount,O.nbnf_nsd_linecount,O.nbnf_etalim_linecount,O.nbnf_ptcut_linecount,\
#                         O.nbnf_all_Nevents,O.nbnf_nsd_Nevents,O.nbnf_etalim_Nevents,O.nbnf_ptcut_Nevents,\
#                         O.nbnf_all,O.nbnf_nsd,O.nbnf_etalim,O.nbnf_ptcut)

#nB_nFs.b_corr_collect_file()
#nB_nFs.nB_nF()
#nB_nFs.nB_nF_plot(nF,nB)
