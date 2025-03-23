### Script written based on Wei-Min Chen's code

library(dplyr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

merged = read.delim("/ceph/CPHG/T1DGC/USERS/dam8mt/data/MEGA/HLA/association_analysis_merged/all/assoc_logistic_chr6_merged.txt", colClasses = c(A1 = "character"))
merged_freq <- read.delim("/ceph/CPHG/T1DGC/USERS/dam8mt/data/MEGA/HLA/association_analysis_merged/all/assoc_logistic_chr6_merged_freq.txt", colClasses = c(A1 = "character", A2 = "character"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate p value when it equals 0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

merged$a = log(2)/log(10) + pnorm(abs(merged$BETA/merged$SE),log.p=TRUE,lower.tail=FALSE)/log(10)
merged$y = floor(merged$a)
merged$x = round(10^(merged$a - merged$y), digits = 3)

merged$P_cal = ifelse(merged$P==0, print(paste(merged$x,sep="","E",merged$y)), print(merged$P))

# with frequency
merged_freq$a = log(2)/log(10) + pnorm(abs(merged_freq$BETA/merged_freq$SE),log.p=TRUE,lower.tail=FALSE)/log(10)
merged_freq$y = floor(merged_freq$a)
merged_freq$x = round(10^(merged_freq$a - merged_freq$y), digits = 3)

merged_freq$P_cal = ifelse(merged_freq$P==0, print(paste(merged_freq$x,sep="","E",merged_freq$y)), print(merged_freq$P))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare final table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

merged = merged %>% 
  select(CHR, SNP, BP, A1, TEST, NMISS, BETA, SE, L95, U95, STAT, P_cal)

colnames(merged)[ncol(merged)] = "P"

# with frequency
merged_freq = merged_freq %>% 
  select(CHR, SNP, BP, A1, A2, TEST, NMISS, BETA, SE, L95, U95, STAT, P_cal, MAF, NCHROBS) %>% 
  rename(P = P_cal)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.table(merged, file = "/ceph/CPHG/T1DGC/USERS/dam8mt/data/MEGA/HLA/association_analysis_merged/all/assoc_logistic_chr6_merged_p.txt", row.names = F, sep = "\t", quote = F)
write.table(merged_freq, file = "/ceph/CPHG/T1DGC/USERS/dam8mt/data/MEGA/HLA/association_analysis_merged/all/assoc_logistic_chr6_merged_freq_p.txt", row.names = F, sep = "\t", quote = F)
