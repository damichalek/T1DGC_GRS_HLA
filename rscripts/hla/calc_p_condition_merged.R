### Script written based on Wei-Min Chen's code

library(dplyr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define arguments
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

args = commandArgs(trailingOnly=TRUE)
val = args[1]
number = args[2]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define variables
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

path = "/ceph/CPHG/T1DGC/USERS/dam8mt/data/MEGA/HLA/association_analysis_merged/conditional_analysis/output/all/"
filename1 = "merged_chr6_"
filename2 = "_condition"
extension = ".assoc.logistic"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

merged = read.table(paste0(path,filename1,val,filename2,number,extension), sep = "", header = T, colClasses = c(A1 = "character"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate p value when it equals 0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

merged$a = log(2)/log(10) + pnorm(abs(merged$BETA/merged$SE),log.p=TRUE,lower.tail=FALSE)/log(10)
merged$y = floor(merged$a)
merged$x = round(10^(merged$a - merged$y), digits = 3)

merged$P_cal = ifelse(merged$P==0, print(paste(merged$x,sep="","E",merged$y)), print(merged$P))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare final table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

merged = merged %>% 
  select(CHR, SNP, BP, A1, TEST, NMISS, BETA, SE, L95, U95, STAT, P_cal)

colnames(merged)[ncol(merged)] = "P"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.table(merged, file = paste0(path,filename1,val,filename2,"_p",number,extension), row.names = F, sep = "\t", quote = F)
