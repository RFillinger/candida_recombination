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
	theme_minimal() + 
	labs( y = "Recombination count", x = "")

ggsave("boxxy_schtufyschtuf.png", width = 7, height = 7 )

