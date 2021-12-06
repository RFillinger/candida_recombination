
df_SC = read.table("aSC5314.SC5314.varAnnot.vcf", sep = "\t")
df_P6 = read.table("aSC5314.P60002.varAnnot.vcf", sep = "\t")

colnames(df_SC) = c("sc.chr","sc.pos","sc.dut1","sc.ref","sc.alt","sc.some.score","sc.pass","sc.info1","sc.info2","sc.info3","sc.info4")
colnames(df_P6) = c("p6.chr","p6.pos","p6.dut1","p6.ref","p6.alt","p6.some.score","p6.pass","p6.info1","p6.info2","p6.info3","p6.info4")

df_SC$"cpos" = paste(df_SC$sc.chr,df_SC$sc.pos)
df_P6$"cpos" = paste(df_P6$p6.chr,df_P6$p6.pos)

aligned_df = merge(df_SC, df_P6, by="cpos")

aligned_df_srt = aligned_df[order(aligned_df[,2],aligned_df[,3]),]

# write.table(aligned_df_srt, file = "sc_p6_shared_SNP_vcf.csv", sep = ",")

df_SC$strain = rep("SC",nrow(df_SC))
df_P6$strain = rep("P6", nrow(df_P6))

colnames(df_P6) = c("chr","pos","dut1","ref","alt","some.score","pass","info1","info2","info3","info4","cpos","strain")
colnames(df_SC) = c("chr","pos","dut1","ref","alt","some.score","pass","info1","info2","info3","info4","cpos","strain")

all = rbind(df_P6, df_SC)

all = all[(all$chr == "Ca21chr1_C_albicans_SC5314")|
		  (all$chr == "Ca21chr2_C_albicans_SC5314")|
		  (all$chr == "Ca21chr3_C_albicans_SC5314")|
		  (all$chr == "Ca21chr4_C_albicans_SC5314")|
		  (all$chr == "Ca21chr5_C_albicans_SC5314")|
		  (all$chr == "Ca21chr6_C_albicans_SC5314")|
		  (all$chr == "Ca21chr7_C_albicans_SC5314")|
		  (all$chr == "Ca21chrR_C_albicans_SC5314"),]

reduced = cbind(all$cpos,all$ref, all$alt, all$info1, all$strain)
colnames(reduced) = c("pos","ref","mutation","info1","strain") # Info1 contains allele frequency

reduced = as.data.frame(reduced)
u_reduced = reduced[!duplicated(reduced$"pos"),] 

# 14,029 unique heterozygous mutations in P60002
# 15,644 unique homozygous mutations in P60002



write.csv(u_reduced, file = "sc_p6_unique_SNP_vcf.csv") 
# Next thing to do is to figure out if all the SNPs in this data frame are unique between SC and P6 (info is in the strain column) It looks off. 