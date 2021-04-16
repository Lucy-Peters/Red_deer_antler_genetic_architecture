

#create list of all window file names 
windows<-list.files(path = "C:/Users/s1767711/Documents/Red_deer/Antler_genetic_architechture/Regional_h2/SNP_windows", pattern = "^window_.*txt$")
windows<-gsub(".txt", "", windows)
windows_list<-as.list(windows)

save(windows_list, file="SNP_windows_list.RData")

load("SNP_windows_list.RData")

x<-readLines("template.R")
f=1

for (l in windows_list){
  write(paste0("w=", "'",l,"'"), paste0("run",f,".R"))
  write(x, paste0("run",f,".R"), append = T)
  f=f+1
}

#make file with list of script names

write(paste0("run1"), paste0("Rscript_list.txt"))

for (i in c (2:3608)){
  write(paste0("run", i), paste0("Rscript_list.txt"),append=T)
}


#make data frame with window ID and corresponding script name
script_list<-read.table("Rscript_list.txt", header=F)
windows<-list.files(path = "C:/Users/s1767711/Documents/Red_deer/Antler_genetic_architechture/Regional_h2/SNP_windows", pattern = "^window_.*txt$")
windows<-gsub(".txt", "", windows)
script_window_list<-data.frame(script.name=script_list$V1, window.name=windows)

write.table(script_window_list, "Rscript_and_window_names.txt", sep="\t", col.names = T, row.names = F,
            quote = F)
