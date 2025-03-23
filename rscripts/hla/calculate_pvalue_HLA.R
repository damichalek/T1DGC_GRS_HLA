### Script written based on Wei-Min Chen's code

library(dplyr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EUR <- read.delim("/ceph/CPHG/T1DGC/USERS/dam8mt/data/MEGA/HLA/association_analysis/Temp/assoc_logistic_chr6_HLA_EUR.txt", colClasses = c(A1 = "character"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate p value when it equals 0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EUR$a = log(2)/log(10) + pnorm(abs(EUR$BETA/EUR$SE),log.p=TRUE,lower.tail=FALSE)/log(10)
EUR$y = floor(EUR$a)
EUR$x = round(10^(EUR$a - EUR$y), digits = 3)

EUR$P_cal = ifelse(EUR$P==0, print(paste(EUR$x,sep="","E",EUR$y)), print(EUR$P))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare final table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EUR = EUR %>% 
  select(CHR, SNP, BP, A1, TEST, NMISS, BETA, SE, L95, U95, STAT, P_cal)

colnames(EUR)[ncol(EUR)] = "P"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.table(EUR, file = "/ceph/CPHG/T1DGC/USERS/dam8mt/data/MEGA/HLA/association_analysis/Temp/assoc_logistic_chr6_HLA_EUR_p.txt", row.names = F, sep = "\t", quote = F)
