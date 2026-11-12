
import pandas as pd
import scrublet as scr
import scipy.io


counts_matrix=pd.read_csv("/home/data/t180330/LAMbrechts_anti_PD1_chemo/0912.write.table.to.srublet/sample21.rdscount.txt",sep="\t",  index_col=0)
counts_matrix=counts_matrix.T
scrub=scr.Scrublet (counts_matrix,expected_doublet_rate=0.05)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
X=counts_matrix
X['doublet']=predicted_doublets
X=X.loc[X['doublet']==False]
X=X.drop(columns=['doublet'])
X.T.to_csv("/home/data/t180330/LAMbrechts_anti_PD1_chemo/0912.write.table.to.srublet/0912.filterf.doubkltr/doubletrilter.sample21.txt",sep="\t")
