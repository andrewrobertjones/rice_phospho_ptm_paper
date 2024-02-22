import os
import pandas as pd
import sys



def read_file(filepath,file):

    locus = file[:-4] #remove csv
    f = open(filepath,"r")
    var_to_stats = {}
    counter = 0
    for line in f:
        line = line[:-1]
        if counter > 0:
            cells = line.split(",")
            variety = cells[0]
            base_counter = 0
            good_snps = 0
            for cell in cells[1:]:
                base_counter += 1
                if cell == "A" or cell == "T" or cell == "G" or cell == "C":
                    good_snps += 1
            df.loc[locus, "locus"] = locus
            #print(variety,good_snps,base_counter)
            df.loc[locus,variety] = good_snps / base_counter
        counter +=1

file_counter = 0

folder = sys.argv[1]
output_file = sys.argv[2]

#folder = "D:/Git/RiceDatabases/Genes-3K-Base-CSV-20211130T093248Z-001/Genes-3K-Base-CSV/MSUgene-Base3K-chr1/1"
for root, dirs, files in os.walk(folder):
    for file in files:
        filepath = os.path.join(root, file)
        #First file, build data frame
        if file_counter == 0:
            all_varieties = ["locus"]
            f_first = open(filepath,"r")
            counter = 0
            for line in f_first:
                line = line[:-1]
                if counter > 0:
                    cells = line.split(",")
                    variety = cells[0]
                    all_varieties.append(variety)
                counter +=1
            df = pd.DataFrame(columns=all_varieties)
            df.set_index("locus")
        file_counter+=1

        if file_counter % 50 == 0:
            print("file count",file_counter)

        read_file(filepath,file)
        #print(df)

df.to_csv(output_file,index=False)
#df.to_csv("D:/Dropbox/DocStore/ProteomicsSoftware/PTMExchange/Rice_build/profile_loci/profiled_loci.csv",index=False)


