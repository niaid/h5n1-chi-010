
##############################################################
# APT command

# create output directory
dir.create(output.dir, showWarnings = FALSE)

apt.path.command=file.path(apt.path,"apt-probeset-summarize")

cel.files=list.files(path=cel.dir,pattern=".cel",full.name=T)
input.command=paste(cel.files,collapse=" ")

affy.anno.pfg=list.files(affy.anno.dir,".pgf",full=TRUE)
affy.anno.mps=list.files(affy.anno.dir,".mps",full=TRUE)
affy.anno.clf=list.files(affy.anno.dir,".clf",full=TRUE)
affy.anno.bgp=list.files(affy.anno.dir,".bgp",full=TRUE)
affy.anno.qcc=list.files(affy.anno.dir,".qcc",full=TRUE)

affy.anno.command=paste("-p",affy.anno.pfg,"-m",affy.anno.mps,"-c",affy.anno.clf,"-b",affy.anno.bgp,"--qc-probesets",affy.anno.qcc,sep=" ")

output.command=paste("-o",output.dir,sep=" ")

system.command=paste(apt.path.command,apt.par,affy.anno.command,output.command,input.command,sep=" ")

# processed cel files by APT
system(system.command)

