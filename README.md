# rice-snp-ptm-analysis
Analysis of SNPs and PTMs for rice

Setup.

Need to download zip file from a) https://drive.google.com/file/d/1ibemgNFD_2tT2rRUfxcRU51ZGnRN12ij and
b) https://drive.google.com/file/d/1i-VKyKe3UqCTDrRUreB3iMqJJp1A6BuF.

These need to be unzipped and referenced following style in run_details.txt.

This reads run_details.txt, which has five parameters:
1. Location of input file (see below for format)
2. Location of where to write output
3. threads=n where n is the number of threads available to the code
4. Location of the SNP-Seek 3K Base SNP download (from above b)
5. Location of the other Rice Databases needed (from above a)

The primary run for the code is done via CreateProteinAlignmentFromSNPSeek.py

The input file should have one RAP-DB or MSU transcript identifier per line, like this:

Os05t0392300-02
Os07t0694000-01
Os01t0549400-01
Os09t0103700-01
Os01t0614500-01
Os01t0588500-01
LOC_Os08g05540.1


The code goes off to SNP-Seek and downloads SNP info in json format. It then extracts the coordinates and bases. It recreates the protein sequence using gff snippets (in Databases) for any given gene, then switches in the base changes to create a new protein molecule.

If the json has been downloaded before, the code uses the cached (local) version.