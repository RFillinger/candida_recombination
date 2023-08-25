library("ggplot2")
library("stringr")
library("multcompView")

main = function(){ 
	# 
	# df = read.csv(file.choose(), header = TRUE, stringsAsFactors = FALSE)

	# df$prop.dead = df$Dead /(df$Alive + df$Dead)
	
	color_1 = "#E7872B" # Orange
	color_2 = "#825CA6" # Purple

	args = commandArgs( trailingOnly = TRUE )

	filename = args[1]

	df = read.csv( filename )
	# print(df)
	df = df[order(df$"Avg..fil..Score"),]

	df$"Image" = factor(df$"Image", levels = df$"Image")

	y = ggplot( data = df, aes(y = as.numeric(`Avg..fil..Score`), x = `Image` )) + 
			geom_bar( stat = "identity", aes(fill = as.factor(`color`))) + 
			theme_bw() +  
			scale_fill_manual( values = c( "#E7872B", "#825CA6", "black" )) + 
			theme(  
					axis.text = element_blank(), 
					axis.title = element_blank(), 
					legend.text = element_blank(),
					axis.text.x = element_blank(),
					panel.grid.major.x = element_blank(),
					panel.grid.minor.x = element_blank(), 
					axis.ticks.x = element_blank(), 
					axis.ticks.y = element_blank(), 
					legend.position = "none",
					panel.border = element_blank() ) + 

			labs( x = "", y = "")

		ggsave( paste0(substr(filename,1,nchar(filename)-4),"_fil_score.png"), height = 7, width = 7 )


}

main()