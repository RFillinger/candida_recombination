library("ggplot2")
library("ggpubr")

dat = read.csv("ddr_wgs_recom_vs_mkr_count.csv")

dat = dat[dat$allele == "n",]

p = ggplot( data = dat, aes( x = `Recombination.count`, y = log(`mkr_count`,2) )) +
	geom_point() + 
	labs(x = "Recombination count", y = "log2(Number of Markers)") 

ggsave("schtufyschtuf.png", width = 7, height = 7 )