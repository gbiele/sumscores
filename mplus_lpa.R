library(data.table)
library(MplusAutomation)
library(tidyLPA)
library(tibble)
source("write_mplus_LPA_model.R")

fnames = c("EF_Q","SO_Q","AL_Q","ER_Q", "LIQ_T","IH_T")


FS = fread("INP/ESEM2018/6F_SR.sav")
FS = FS[,tail(1:(ncol(FS)),length(fnames)), with =  F]
colnames(FS) = fnames
FS = apply(FS,2,scale)
load("INP/ESEM2018/full_my_data.Rdata")
my_data[,ADHD := factor(ADHD, labels = c("No","Yes"))]

FSadhd = FS[my_data$ADHD == "Yes",]
FSadhd = FS

prepareMplusData(data.frame(FSadhd),"LPA/dd.dat")

bc = ""
for (m in c(2,4,5,6)) {
  for (p in 2:5) {
    write_mplus_LPA_model(as.tibble(FSadhd),
                          EF_Q, SO_Q, AL_Q, ER_Q, LIQ_T, IH_T,
                          n_profiles = p,
                          model = m,
                          script_filename = paste0("LPA/m",m,"_p",p,".inp"),
                          output_filename = paste0("LPA/m",m,"_p",p,".out"),
                          savedata_filename = paste0("m",m,"_p",p,"_sd.out"),
                          data_filename = "LPA/dd.dat",
                          starts = c(4000*m, 400*m), m_iterations = 3000*m, st_iterations = 20*m,
                          n_processors = 4,
                          remove_tmp_files = F)
    
    bc = paste0(bc,"\n mplus m",m,"_p",p,".inp")
  }
}

cat(bc,file = "LPA/batch.bat")