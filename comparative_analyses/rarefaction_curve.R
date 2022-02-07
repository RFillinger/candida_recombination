library("ggplot2")
library("ggpubr")

dat = read.csv("ddr_wgs_recom_vs_mkr_count.csv")

dat = dat[dat$allele == "n",]

p = ggplot( data = dat, aes( x = log(`mkr_count`,2), y = `Recombination.count`  )) +
	geom_point() + 
	labs(y = "Recombination count", x = "log2(Number of Markers)") + 
	theme_minimal() + 
	stat_smooth( method = "lm", formula = y ~ poly(x,2), size = 1, se = FALSE, color = "red" )

ggsave("schtufyschtuf.png", width = 7, height = 7 )


w = ggplot( data = dat, aes( x = `data.source`, y = `Recombination.count`, group = `data.source`)) + 
	geom_boxplot( fill = "skyblue", notch = FALSE ) + 
	geom_jitter( size = 1, color = "black", width = 0.3 ) + 
	theme_bw() + 
	ylim(c(0,600)) + 
	# labs( y = "Recombination count", x = "")
	labs( y = "", x = "") + 
	theme( panel.grid.major.x = element_blank(), 
		panel.grid.minor.x = element_blank(), 
		panel.grid.minor.y = element_blank(), 
		panel.border = element_blank(),
		axis.title.x = element_blank(), 
		axis.ticks.x = element_blank(), 
		axis.text.x = element_blank(),
		axis.title.y = element_blank(), 
		axis.ticks.y = element_blank(), 
		axis.text.y = element_blank() ) 

ggsave("boxxy_schtufyschtuf.png", width = 7, height = 7 )

