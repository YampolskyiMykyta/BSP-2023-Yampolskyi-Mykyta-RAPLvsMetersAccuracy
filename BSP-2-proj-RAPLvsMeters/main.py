import readdata2 as rd2

rd2.rewriteFileMETERandDeltas("dataMETER.csv","dataDELTAS.csv",False)
rd2.rewriteFileMETERandDeltas("dataMETER_ver2.csv", "dataDELTAS_ver2.csv", True)  # True - less conservative approach
#rd2.get_everything_about_RMD()

rd2.get_speedup_rapl_median()
#rd2.show_plot(True)

#rd2.get_data_before_and_after()

