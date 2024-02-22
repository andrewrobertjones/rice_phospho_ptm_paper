from os import listdir
from os.path import join
import os
import shutil

#path_to_SNP_files = "D:/Git/RiceDatabases/remaining-genes-csv/remaining-genes-csv/"
path_to_SNP_files = "/home/arjones/hc-storage/DataSets/RicePhosphoSNPs/remaining-genes-csv/"
#output_folder = "D:/Git/RiceDatabases/Genes-3K-Base-CSV-20211130T093248Z-001/Genes-3K-Base-CSV/MSUgene-Base3K-chr1/"
output_folder = "/home/arjones/hc-storage/DataSets/RicePhosphoSNPs/Genes-3K-Base-CSV/"
#1/LOC_Os01g01390/"

for f in listdir(path_to_SNP_files):

    filename = join(path_to_SNP_files, f)
    if ".csv" in filename and "LOC_" in filename:
        #print("Chr:",f[6:8],"chr_int:",)
        chr_num = int(f[6:8])   #LOC_Os01g38510.1.csv grab 01 and covert to 1

        locus_path = output_folder + str(chr_num) + "/" + f.split(".")[0] + "/"
        if not os.path.exists(locus_path):
            os.makedirs(locus_path)

        newfilename =  locus_path + f
        shutil.copyfile(filename, newfilename)
        print("writing to ",newfilename)







