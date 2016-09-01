import qgsm_collect as collect
import qgsm_FB_corr as fb_corr

energy = 900 
#O = collect.QGSM_Distributions(energy)
#O.collectData()
#O.writeAnalysis()

FBs = fb_corr.FB(energy)


FBs.b_corr_collect_file()
#FBs.b_corr_collect_mem(O.ALL,O.lineCount,O.countedEvents)
FBs.b_corr_count()
FBs.b_corr_plot()
