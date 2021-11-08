
library("ggplot2")
library("stringr")

# Command: 
# Rscript old_program_versions/cf_graph.R 

path_components <- function( path_and_file ){

	file_path_vector = unlist(str_split( path_and_file, "/" ))
	filename = file_path_vector[length(file_path_vector)]
	path_components = file_path_vector[1:length(file_path_vector)-1] 
	path_string = paste0( path_components, "/" , collapse = "/" )

	file_stuff <- c( path_string, filename )

	return( file_stuff )

}

main <- function() {

	args = commandArgs( trailingOnly = TRUE )

	file_vec = c( "SCxP60002/allele_files_for_ploidy_measurements/all_inclusive_markers.csv", 
				"SCxP60002/allele_files_for_ploidy_measurements/excluded_markers.csv", 
				"SCx529L/allele_files_for_ploidy_measurements/all_inclusive_markers.csv", 
				"SCx529L/allele_files_for_ploidy_measurements/excluded_markers.csv" )

	for ( filename in file_vec ){

		df = read.csv( filename, header = FALSE )
		colnames( df ) = c( "sample.name", "cat.num", "chromosome", "chr.loc", "allele.count", "mkr.depth", "allele.1.prop", "allele.2.prop")
		df$"allele.count" = as.factor( df$"allele.count")
		df = df[ as.numeric(df$"allele.count") >= 2, ]
		df_cutoff = df[ (df$"allele.1.prop" < 55)&(df$"allele.1.prop" > 45)&(df$"allele.2.prop" > 45)&(df$"allele.2.prop" < 55), ]

		color_vec_length = as.numeric(length(unique( df$"allele.count" )))
		# color_vec = c( "black", "red", "blue", "green", rep("black", color_vec_length-4))
		color_vec = c( "red", "blue", "green", rep("black", color_vec_length-3))

		print( paste0( "Total number of heterozygous markers: ", nrow(df[(df$"mkr.depth"),])))
		print( paste0( "Number of usable heterozygous markers: ", nrow(df_cutoff)))
		
		p = ggplot( df, aes( x= `allele.1.prop`, y = `allele.2.prop`, colour = `allele.count`) ) + 
			lims(x = c(0,100), y = c(0,100)) + 
			geom_point() + 
			theme_minimal() + 
			scale_colour_manual( values = color_vec ) + 
			labs( x = "Major allele proportion", y = "Minor allele proportion ", colour = "Number of alleles in marker")

		ggsave( paste0( substr(filename, 1, str_length(filename)-4), ".png" ), width = 7, height = 7 )

		c = ggplot( df, aes( x = `allele.1.prop`, y = `allele.2.prop` ))+ 
			geom_density_2d() +
			lims(x = c(0,100), y = c(0,100)) + 
			theme_minimal() + 
			labs( x = "Major allele proportion", y = "Minor allele proportion ", colour = "Number of alleles in marker")
			

		ggsave( paste0( substr(filename, 1, str_length(filename)-4), ".contour.png" ), width = 7, height = 7 )

		# g = ggplot( df_cutoff, aes( x= `allele.1.prop`, y = `allele.2.prop`, colour = `allele.count`) ) + 
		# 	lims(x = c(0,100), y = c(0,55)) +
		# 	geom_point() + 
		# 	theme_minimal() + 
		# 	scale_colour_manual( values = color_vec ) +
		# 	labs( x = "Major allele proportion", y = "Minor allele proportion", colour = "Number of alleles in marker")

		# ggsave( paste0( substr(filename, 1, str_length(filename)-4), ".cutoff.png" ), width = 7, height = 7 )

	}

}

main()