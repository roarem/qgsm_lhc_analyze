import qgsm_collect as collect
import qgsm_FB_corr as fb_corr

energy = 900
O = collect.QGSM_Distributions(energy)
O.collectData(bcorr=False,nbnf=True)
O.create_hists()
#O.writeAnalysis()

#FBs = fb_corr.FB(energy)
#FBs.b_corr_collect_file()
#FBs.b_corr_collect_mem(O.bcorr_Nevents,O.bcorr,O.bcorr_linecount)
#FBs.b_corr_count(O.bcorr_Nevents,O.bcorr,O.bcorr_linecount)
#FBs.b_corr_plot(from_file=False)

#nB_nFs = fb_corr.FB(energy)
#nB_nFs.nB_nF_collect_mem(O.nbnf_all_linecount,O.nbnf_nsd_linecount,O.nbnf_etalim_linecount,O.nbnf_ptcut_linecount,\
#                         O.nbnf_all_Nevents,O.nbnf_nsd_Nevents,O.nbnf_etalim_Nevents,O.nbnf_ptcut_Nevents,\
#                         O.nbnf_all,O.nbnf_nsd,O.nbnf_etalim,O.nbnf_ptcut)

#nB_nFs.b_corr_collect_file()
#nB_nFs.nB_nF()
#nB_nFs.nB_nF_plot(nF,nB)
