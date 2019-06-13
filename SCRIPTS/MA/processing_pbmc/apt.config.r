# this module cannot work under interactive node
print("This module does not work under interactive node")

##############################################################
# input parameters

# path for cel files
cel.dir="/sg_data/CHI/PROJECTS/H5N1/MICROARRAY/LEVEL1/2014_11_24/HuGene-2_1-st/CELfiles/"

# output path
output.dir="./gexp_APT"

# set up path for APT
apt.path="/usr/local/bio_apps/apt/bin/"

# rma option
apt.par="-a rma-sketch"

# pdg, mps, clf, bgp and qcc files
affy.anno.dir="/sg_data/PROJECTS/ARRAY_LIB/HuGene-2_1-st"


