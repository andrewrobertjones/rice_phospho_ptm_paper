
f = open("D:/Dropbox/DocStore/Articles/2022_rice_phosphoproteome/SNP_analysis/MSU_SNP_all.txt","r")

counter = 1

output_folder  = "D:/Dropbox/DocStore/Articles/2022_rice_phosphoproteome/SNP_analysis/slurm_run_details/"

text_for_run_details = ""

location_for_run = "/home/arjones/hc-storage/DataSets/RicePhosphoSNPs/rice-snp/rice_snp_analysis/"
text_for_run_details = "/home/arjones/hc-storage/DataSets/RicePhosphoSNPs/rice-snp/output_SAAVs/\n"
text_for_run_details += "threads=20\n"
text_for_run_details += "/home/arjones/hc-storage/DataSets/RicePhosphoSNPs/Genes-3K-Base-CSV/\n"
text_for_run_details += "/home/arjones/hc-storage/DataSets/RicePhosphoSNPs/\n"

current_transcripts = ""
for line in f:
    transcript_id = line.rstrip()
    current_transcripts += transcript_id + "\n"
    print("c:",counter)
    if counter % 50 == 0:
        print("here")
        snp_file_name = "MSU_SNPs"+str(counter)+".txt"
        f_out2 = open(output_folder+snp_file_name,"w")
        f_out2.write(current_transcripts)
        f_out2.close()
        current_transcripts = ""

        f_out = open(output_folder+"run_details"+str(counter)+".txt","w")
        f_out.write(location_for_run + snp_file_name+"\n")
        f_out.write(text_for_run_details)
        f_out.close()
    counter+=1

#Last one
snp_file_name = "MSU_SNPs" + str(counter) + ".txt"
f_out2 = open(snp_file_name, "w")
f_out2.write(current_transcripts)
f_out2.close()

f_out = open("run_details" + str(counter) + ".txt", "w")
f_out.write(location_for_run + "snp_file_name\n")
f_out.write(text_for_run_details)
f_out.close()




