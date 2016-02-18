##merge2files

Script to compare and merge two files based on user specified matched columns.
for usage:
merge2files.pl -h

##runVirusdetection
* Tool to detect virus or bacterial sequences in sequenced samples using GOTTCHA
* Tracey Allen K. Freitas, Po-E Li, Matthew B. Scholz and Patrick S. G. Chain (2015) Accurate read-based metagenome characterization using a hierarchical suite of unique signatures, Nucleic Acids Research (DOI: 10.1093/nar/gkv180)
* Download latest signature database for specific organism from : ftp://ftp.lanl.gov/public/genome/gottcha/GOTTCHA_database_v20150825/
* Also needs lookup database GOTTCHA_lookup.tar.gz
* runVirusdetection.sh 
##runPathogenDetection 
* This is generic perl wrapper to detect pathogen seuences in using GOTTCHA method [see above description for reference]
* More robust options to select unmapped reads sequences
