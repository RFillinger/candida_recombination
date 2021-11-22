library("car")
library("dgof")
# Command to run this program: Rscript recom_stats.R SCx529L/recom_529LxSC_4way.csv SCxP60002/recom_P60002xSC_4way.csv

main <- function() {

	args = commandArgs( trailingOnly = TRUE )
	recom_stats_file_1 = args[1]
	recom_stats_file_2 = args[2]

	df_1 = read.csv( recom_stats_file_1 )
	df_2 = read.csv( recom_stats_file_2 )

	# Combine the two crosses into a single column of a data frame, and create a separate column denoting which value belongs to which cross. 
	recom.ratio = c(df_1$recom.by.mkr, df_2$recom.by.mkr)
	group = c( rep("1", length(df_1$recom.by.mkr)), rep("2", length(df_2$recom.by.mkr)) )
	
	df = data.frame(recom.ratio, group)

	leveneTest( df$recom.ratio, df$group ) # If p-value is significant, it indicates that the two groups ahve different variances 

	ks.test( df_1$recom.by.mkr, df_2$recom.by.mkr ) 

}

main()
