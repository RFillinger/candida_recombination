
library("stringr")
library("ggplot2")
library("grid")
library("chisq.posthoc.test")
library("ggpubr")

file_extension = ".png"

path_components <- function( path_and_file ){

	file_path_vector = unlist(str_split( path_and_file, "/" ))
	filename = file_path_vector[length(file_path_vector)]
	path_components = file_path_vector[1:length(file_path_vector)-1] 
	path_string = paste0( path_components, "/" , collapse = "/" )

	file_stuff <- c( path_string, filename )

	return( file_stuff )

}


absolute_positions <- function( chr_lengths ){

	univ_pos_names = names( chr_lengths )
	univ_pos = vector( mode="list", length=length( chr_lengths ) )
	total = 0 
	i = 1
	for ( chromosomes in chr_lengths ){ 
		univ_pos[[ i ]] <- total 
		i = i + 1
		total = total + chromosomes 
	}
	names( univ_pos ) <- univ_pos_names 
	return( univ_pos )

}


geno_graphs <- function( path_and_file, chr_labs, chr_lengths, universal_positions,
	ymax = 250, ignore_missing = 1, tick_perct = 0.025, strangers = 0,
	color_1 = "#E7872B", color_2 = "#825CA6", color_3 = "#5AAA46", theme_color = "black" ) {

	path = path_components( path_and_file )[1]
	file = path_components( path_and_file )[2]

	df = read.csv( path_and_file )

	if (strangers){
		allele_sets <- c("a", "b", "h", "missing", "s")
	} else {
		allele_sets <- c("a", "b", "h", "missing")
	}

	chromosomes <- unique( df$chr )
	max_sizes = vector( mode="list", length=length( chromosomes ) )
	names( max_sizes ) = as.vector(chromosomes)
	for ( chromosome in chromosomes ){
		# Subset the data frame to obtain the maximum sizes of chromosomes in the df
		chr_subset = df[df$chr == chromosome, ]
		max_sizes[[chromosome]] <- max(chr_subset$univ.pos) 
	}

	chromosomes_delineations = unlist(max_sizes)[c(1:length(chromosomes)-1)] 
	for ( allele in allele_sets ){

		if ( allele == "a" ) {
			y_var = df$a
			color = color_1
		} else if ( allele == "b" ) {
			y_var = df$b
			color = color_2
		} else if ( allele == "h" ){
			y_var = df$h
			color = color_3
		} else if ( allele == "missing" ){
			y_var = df$missing
			color = "gray"
		} else if ( (allele == "s") & strangers ){
			y_var = df$s
			color = "black"
		}
		
		# Make the plots 
		p <- ggplot( data = df, aes( x = `univ.pos`, y = y_var )) + 
		geom_col( color = color ) + 
		# geom_vline( xintercept = chromosomes_delineations, size = 1, color = "black" ) +
		theme_minimal() + 
		theme( legend.position = "none", axis.text.x = element_blank() ) + 
		labs( y = "Allele Proportion", x = NULL ) + 
		ylim( c(0,ymax))	

		# One loop for the delineations 
		for ( chromosomes in names(chromosomes_delineations) ){
			x_spot = chr_lengths[[chromosomes]] + universal_positions[[ chromosomes ]]
			p <- p + geom_segment( x = x_spot, y = 0, xend = x_spot, yend = -tick_perct * ymax, color = theme_color)
		}

		# Second loop for the labels that go under the delineations
		i = 1 
		for ( chromosomes in names(chr_labs) ){
			x_spot = round(chr_lengths[[ chromosomes ]]/2 + universal_positions[[ chromosomes ]])
			p <- p + annotation_custom( textGrob(chr_labs[[chromosomes]], gp = gpar( fontsize = 10, col = theme_color)), # Adds the chromosome labels
							xmin = x_spot , xmax = x_spot, ymin = -tick_perct/2 * ymax, ymax = -tick_perct/2 * ymax )
			i = i + 1
		}

		ggsave( paste0(path, allele, "_", substr( "genotype_tally", 1, str_length(file)-4 ), file_extension ), width=7, height=7 )

	}

}


recom_graphs <- function( path_and_file, chr_labs, chr_lengths, universal_positions,  
							ymax = 1.5, tick_perct = 0.025, add_stdev_lines = 0, color = "gold", theme_color = "black" ) {

	path = path_components( path_and_file )[1]
	file = path_components( path_and_file )[2]

	df = read.csv( path_and_file )

	# Chi-squared test (goodness of fit) lines go here
	equal_probability_vector = rep((1/length(df$recom.by.mkr)), length(df$recom.by.mkr)) # Assumes random chance of a recombination for each bin. 
	redistribution_vector = rep( round(sum(df$recom.by.mkr)/length(df$recom.by.mkr), 0), length(df$recom.by.mkr) )

	chi_squared = chisq.test( x = df$recom.by.mkr, p = equal_probability_vector )
	# print( chi_squared$expected  ) # Every value should be above 5

	shap_test = shapiro.test(df$recom.by.mkr) # Test for normalcy of distributions

	recom_stdev = sd(df$recom.by.mkr)
	recom_avg = mean(df$recom.by.mkr)
	threshold = recom_avg + (2*recom_stdev) # represents 95% of data
	opp_threshold = -threshold

	big_bins = df[(df$recom.by.mkr>threshold)|(df$recom.by.mkr<opp_threshold),]

	y_pos = -0.5
	col_1 = "magenta"
	col_2 = "turquoise"

	histogram <- ggplot(df, aes(x=df$recom.by.mkr)) + 
		geom_histogram(binwidth=1) + 
		lims( x=c(0,150) ) +
		labs( y = "Frequency", x = "Recombination events") + 
		geom_segment( x = recom_avg, y = 0, xend = recom_avg, yend = y_pos, color = col_1) + 
		geom_segment( x = threshold, y = 0, xend = threshold, yend = y_pos, color = col_2 ) + 
		geom_segment( x = opp_threshold, y = 0, xend = opp_threshold, yend = y_pos, color = col_2 ) + 
		annotation_custom( textGrob("Mean", gp = gpar( fontsize = 10, col = col_1)), # Adds label to mean line
							xmin = recom_avg , xmax = recom_avg, ymin = y_pos, ymax = y_pos ) + 
		annotation_custom( textGrob("Mean + 2σ", gp = gpar( fontsize = 10, col = col_2 )), # Adds label to mean line
							xmin = threshold , xmax = threshold, ymin = y_pos, ymax = y_pos )


	ggsave( paste0(path, substr( file, 1, str_length(file)-4 ), ".recom_hist", file_extension ), width=7, height=7 )

	chromosomes <- unique( df$chr )
	max_sizes = vector( mode="list", length=length( chromosomes ) )
	names( max_sizes ) = as.vector(chromosomes)
	for ( chromosome in chromosomes ){
		# Subset the data frame to obtain the maximum sizes of chromosomes in the df
		chr_subset = df[df$chr == chromosome, ]
		max_sizes[[chromosome]] <- max(chr_subset$univ.pos) 
	}

	chromosomes_delineations = unlist(max_sizes)[c(1:length(chromosomes)-1)] 

	y_disp = -(ymax*tick_perct) + 15

	p <- ggplot( data = df, aes( x = `univ.pos`, y = `recom.by.mkr` )  ) + 
		geom_col( color = "cyan" ) + 
		# geom_vline( xintercept = chromosomes_delineations, size = 0.2, color = color ) +
		theme_minimal() + 
		theme( legend.position = "none", axis.text.x = element_blank() ) + 
		labs( y = "Recombinations per Marker", x = NULL ) + 
		ylim( c(0,ymax)) 

	if ( add_stdev_lines ){
		p <- p + geom_hline( yintercept = recom_avg ) + 
				 geom_hline( yintercept = threshold )

	# 	p <- p + 
	# 	geom_segment( x = min(unlist(max_sizes)), y = recom_avg, xend = max(unlist(max_sizes)), yend = recom_avg, color = col_1) + 
	# 	geom_segment( x = min(unlist(max_sizes)), y = threshold, xend = max(unlist(max_sizes)), yend = threshold, color = col_2) + 
	# 	annotation_custom( textGrob("Mean", gp = gpar( fontsize = 10, col = col_1)), # Adds label to mean line
	# 						xmin = 0 , xmax = 0, ymin = recom_avg + (ymax*tick_perct), ymax = recom_avg + (ymax*tick_perct) ) + 
	# 	annotation_custom( textGrob("Mean + 2σ", gp = gpar( fontsize = 10, col = col_2 )), # Adds label to mean line
	# 						xmin = 0 , xmax = 0, ymin = threshold + (ymax*tick_perct), ymax = threshold + (ymax*tick_perct) )
	}

	# One loop for the delineations 
	for ( chromosomes in names(chromosomes_delineations) ){
		x_spot = chr_lengths[[chromosomes]] + universal_positions[[ chromosomes ]]
		p <- p + geom_segment( x = x_spot, y = 0, xend = x_spot, yend = -tick_perct * ymax, color = theme_color)
	}
	

	# Second loop for the labels that go under the delineations
	i = 1 
	for ( chromosomes in names(chr_labs) ){
		x_spot = round(chr_lengths[[ chromosomes ]]/2 + universal_positions[[ chromosomes ]])
		p <- p + annotation_custom( textGrob(chr_labs[[chromosomes]], gp = gpar( fontsize = 10, col = theme_color)), # Adds the chromosome labels
						xmin = x_spot , xmax = x_spot, ymin = -tick_perct * ymax, ymax = -tick_perct * ymax )
		i = i + 1
	}

	if ( add_stdev_lines ){
		ggsave( paste0(path, substr( file, 1, str_length(file)-4 ), "_STDEV", file_extension), width=7, height=7 )
	} else {
		ggsave( paste0(path, substr( file, 1, str_length(file)-4 ), file_extension), width=7, height=7 )
	}
	

}


modelling_data_prep <- function( path_and_file, chr_lengths ) {

	# Open the file and read in the data to a dataframe 
	path = path_components( path_and_file )[1]
	file = path_components( path_and_file )[2]
	df = read.csv( path_and_file )
	df = df[ order( df[,3]), ]

	arm_names = unique( df$chromosome_arm )

	chromosome_arms <- vector( mode="list", length=length( chr_lengths )*2 )

	index = 1
	new_arm_names = c()

	for (chromosomes in names(chr_lengths) ){

		for (arms in arm_names) {
			# Create arm names to store the data in a mapped list and store them in a vector
			arm_name = paste0( chromosomes, "_", arms )
			new_arm_names = c( new_arm_names, arm_name )

			# Find the subset of data that matches both left and right arms of each chromosome
			chr_subset = df[ (df$chromosome_arm == arms), ] # Matches the arm (default is left or right)
			chr_arm_subset = chr_subset[ (chr_subset$chromosome == chromosomes), ] # Matches the chromosome and store them as a subset
			chromosome_arms[[ index ]] = chr_arm_subset # Put it in the list 
			index = index + 1
		}
	}

	# Combine the names with the cumulative recombination data create the final list structure. 
	names( chromosome_arms ) = new_arm_names 
	return ( chromosome_arms )

}


lmp <- function (modelobject) {

    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    if ( !is.null(f) ){
    	p <- pf(f[1],f[2],f[3],lower.tail=F) # This line will error if there are no recombinations in a given chromosome 
    } else {
		p <- "NULL"
	}
    attributes(p) <- NULL
    return(p)

}


centro_graphs <- function( path_and_file, centromeres, chr_lengths, universal_positions, chr_labs, model_dat, 
							ymaximum = 200, tick_perct = 0.01, theme_color = "black" ) {

	path = path_components( path_and_file )[1]
	file = path_components( path_and_file )[2]
	
	df = read.csv( path_and_file )

	df = df[ order( df[,3]), ]
	
	centro_means <- vector( mode="list", length=length( centromeres ))
	names( centro_means ) = as.vector( names(centromeres) )
	for ( centromere in names(centromeres) ){
		mean = round(mean( centromeres[[centromere]] ), digits = 0)
		centro_means[[ centromere ]] <- mean + universal_positions[[centromere]]
	}

	chromosomes <- unique( df$chromosome )
	max_sizes = vector( mode="list", length=length( chromosomes ) )
	names( max_sizes ) = as.vector(chromosomes)
	for ( chromosome in chromosomes ){
		# Subset the data frame to obtain the maximum sizes of chromosomes in the df
		chr_subset = df[df$chromosome == chromosome, ]
		max_sizes[[chromosome]] <- max(chr_subset$absolute_pos) 
	}

	p <- ggplot( data = df, aes( x = `absolute_pos`, y = `cumulative_recom_cnt` )  ) + 
		geom_step( aes(color = `chromosome`) ) + 
		theme_minimal( ) +  
		theme( legend.position = "none", axis.text.x = element_blank() ) + 
		ylim( c(0,ymaximum)) + 
		labs( x = NULL, y = "Cumulative Recombination" )
		
	# I hate that I have to do this... it adds ticks and lines and editable labels to the graph where the centromeres are. Ugh.
	y_disp = -(ymaximum*tick_perct) - 5 

	i = 1 
	for ( centromere in centro_means ){
		p <- p + geom_segment( x = centromere, y = 0, xend = centromere, yend = tick_perct * ymaximum, color = theme_color ) + # Adds the better custom ticks
				annotation_custom( textGrob(chr_labs[[i]], gp = gpar( fontsize = 10, col = theme_color)), # Adds the chromosome labels
					xmin = centromere , xmax = centromere, ymin = y_disp, ymax = y_disp ) 
				# annotation_custom(segmentsGrob(gp = gpar(col = theme_color, lwd = 2)), # Adds the custom ticks
        			# xmin = centromere, xmax = centromere, ymin = y_disp, ymax = y_disp)
		i = i + 1
	}

	ggsave( paste0(path, substr( file, 1, str_length(file)-4 ), file_extension), width=7, height=7 )

	i = 0
	for (arms in names(model_dat)) {

		chr_name = unique( model_dat[[ arms ]]$chromosome )

		if (i%%2 == 1){
			x_loc = round(chr_lengths[[chr_name]]*(1/4)) + universal_positions[[chr_name]]
			y_loc = ymaximum*0.9

		} else {
			x_loc = round(chr_lengths[[chr_name]]*(3/4)) + universal_positions[[chr_name]]
			y_loc = ymaximum*0.7
		}

		arm_model = lm (model_dat[[ arms ]]$cumulative_recom_cnt ~ model_dat[[ arms ]]$absolute_pos)
		r2_val = round(summary(arm_model)$r.squared, 2 )
		model_p_val = lmp( arm_model )
		coeff_p_val = round( summary( arm_model)$coefficients[,4], 2)

		lm_model = model_dat[[arms]][,4] ~ model_dat[[arms]][,3]
		
		# Prints p-values of model fit to each arm
		# print( arms )
		# print( model_p_val )
		
		# This set will generate lines, but make a lot of warnings
		# p = p + geom_smooth( data = model_dat[[arms]], method = "lm", se = FALSE, color = "black", na.rm = TRUE ) + 
		# 		annotation_custom( textGrob( r2_val, gp = gpar( fontsize = 7, col = theme_color)), # Adds the r squared values
		# 			xmin = x_loc , xmax = x_loc, ymin = y_loc, ymax = y_loc )

		# This chunk won't make lines or warnings
		p = p + geom_smooth( data = model_dat[[arms]], formula = lm_model, method = "lm", se = FALSE, color = "black", na.rm = TRUE ) + 
				annotation_custom( textGrob( r2_val, gp = gpar( fontsize = 7, col = theme_color)), # Adds the r squared values
					xmin = x_loc , xmax = x_loc, ymin = y_loc, ymax = y_loc ) #+ 
				# annotation_custom( textGrob( model_p_val, gp = gpar( fontsize = 7, col = theme_color)), # Adds the p-values
				# 	xmin = x_loc, xmax = x_loc, ymin = y_loc - 0.05*ymaximum, ymax = y_loc - 0.05*ymaximum )
		i = i + 1
	}

	ggsave( paste0(path, substr( file, 1, str_length(file)-4 ), "_lm", file_extension ), width=7, height=7 )

}


shared_geno_graphs <- function( path_and_file, chr_labs, chr_lengths, univ_pos, a = "1", b = "2", h = "n", missing = "-", include_missing = 0, 
								color_1 = "#E7872B", color_2 = "#825CA6", color_3 = "#5AAA46" ){

	path = path_components( path_and_file )[1]
	file = path_components( path_and_file )[2]
	df = read.csv( path_and_file )

	for ( leels in unique(df[,2]) ){

		if ( leels == a ) {
			# Sort a first, then b, then h
			df = df[ order(df[,4]), ]
			df[,2] = factor( unique(df[,2]), levels = c(b,h,a) )
			colors = c(color_2, color_3, color_1)

		} else if ( leels == b ) {

			df = df[ order(df[,5]), ]
			df[,2] = factor( unique(df[,2]), levels = c(a,h,b) )
			colors = c(color_1, color_3, color_2)			

		} else if ( leels == h ){

			df = df[ order(df[,6]), ]
			df[,2] = factor( unique(df[,2]), levels = c(a,b,h) )
			colors = c(color_1, color_2, color_3)

		} else if ( (leels == missing) && include_missing ){

			df = df[ order(df[,7]), ]
			df[,2] = factor( unique(df[,2]), levels = c(a,h,b,missing) )
			colors = c(color_1, color_3, color_2, "magenta")

		}

		p <- ggplot( data = df, aes( fill = `allele`, x = `Progeny.Strain`, y = `count` )  ) +
			geom_bar( position = "fill", stat = "identity", width = 2.5 ) + 
			scale_fill_manual( name = "Allele", values = colors, labels = c( "Parent 2", "Heterozygous", "Parent 1") ) + 
			scale_x_discrete( limits = df$Progeny.Strain ) + # Adds parent and progeny names 
			theme_minimal() + 
			labs( x = NULL, y = "Allele Proportion", color = "Allele" ) +
			theme(  axis.text=element_text(size=12), 
					axis.title = element_text(size=14), 
					legend.text = element_text(size = 12),
					axis.text.x = element_blank() )
					# axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0) ) + 
			

		ggsave( paste0(path, substr( file, 1, str_length(file)-4 ), "_", leels , file_extension ), width=10, height=7)

	}
}


track_lengths <- function ( path_and_file, chr_labs, chr_lengths, univ_pos, bandwidth = 25, ymaximum = 600, xmaximum = 2000 ){

	path = path_components( path_and_file )[1]
	file = path_components( path_and_file )[2]
	df = read.csv( path_and_file )

	g <- ggplot( data = df, aes( x = df[,4]/1000)) + 
	geom_histogram( binwidth = bandwidth, color = "black", fill = "white") +
	labs( x = paste0("Estimated Recombination Size (Kbp) Binsize: ", bandwidth), y = "Frequency") +
	ylim( c(0,ymaximum)) +
	xlim( c(0,xmaximum)) + 
	theme_minimal()

	ggsave( paste0( path, substr( file, 1, str_length(file)-4 ), file_extension ), width=7, height=7 )

}

recom_histo <- function( rhf_filename, path_and_file, bandwidth ) {

	path_stuff = path_components( path_and_file )
	path = path_stuff[1]
	file = path_stuff[2]

	df = read.csv( rhf_filename ) 
	x <- ggplot( data = df, aes( x = df[,2] )) + 
	geom_histogram( binwidth = bandwidth, color = "black", fill = "white" ) + 
	labs( x = "", y = "Recombination Frequency" ) + 
	theme_minimal()

	ggsave( paste0( path, "RECO_histo_", substr(file,12,str_length(file)-4 ), ".png" ), width=7, height=7 )

}


main <- function() {

	args = commandArgs( trailingOnly = TRUE )
	chromosome_filename = args[1]
	path_and_file = args[2]
	graphs = args[3]

	chr_df = read.csv( chromosome_filename )

	chromosomes = unique( chr_df$chr_name )
	chr_labs    = vector( mode="list", length=length( chromosomes ))
	chr_lengths = vector( mode="list", length=length( chromosomes ))
	centromeres = vector( mode="list", length=length( chromosomes ))

	names( chr_labs )    = as.vector(chromosomes)
	names( chr_lengths ) = as.vector(chromosomes)
	names( centromeres ) = as.vector(chromosomes)

	for ( chromosome in chromosomes ){
		# Subset the data frame to obtain the maximum sizes of chromosomes in the df
		temp_df = chr_df[chr_df$chr_name == chromosome, ]
		
		chr_labs[[chromosome]]    = temp_df$alias
		chr_lengths[[chromosome]] = temp_df$chr_length
		centromeres[[chromosome]] = c( temp_df$centromere_1, temp_df$centromere_2 )
	}

	# chr_labs = list("Ca21chr1_C_albicans_SC5314" = "1", 
	# 				"Ca21chr2_C_albicans_SC5314" = "2",
	# 				"Ca21chr3_C_albicans_SC5314" = "3",
	# 				"Ca21chr4_C_albicans_SC5314" = "4",
	# 				"Ca21chr5_C_albicans_SC5314" = "5",
	# 				"Ca21chr6_C_albicans_SC5314" = "6",
	# 				"Ca21chr7_C_albicans_SC5314" = "7",
	# 				"Ca21chrR_C_albicans_SC5314" = "R")

	# chr_lengths = list( "Ca21chr1_C_albicans_SC5314" = 3190000,
	# 					"Ca21chr2_C_albicans_SC5314" = 2230000,
	# 					"Ca21chr3_C_albicans_SC5314" = 1800000,
	# 					"Ca21chr4_C_albicans_SC5314" = 1600000,
	# 					"Ca21chr5_C_albicans_SC5314" = 1190000,
	# 					"Ca21chr6_C_albicans_SC5314" = 1030000,
	# 					"Ca21chr7_C_albicans_SC5314" = 950000,
	# 					"Ca21chrR_C_albicans_SC5314" = 2290000 )

	# # Centromere names need to be identical to chromosome names. 
	# centromeres = list( "Ca21chr1_C_albicans_SC5314" = c(1563038,1565967), 
	# 					"Ca21chr2_C_albicans_SC5314" = c(1927255,1930214),
	# 					"Ca21chr3_C_albicans_SC5314" = c(823333,826481),
	# 					"Ca21chr4_C_albicans_SC5314" = c(992579,996216),
	# 					"Ca21chr5_C_albicans_SC5314" = c(468716,471745),
	# 					"Ca21chr6_C_albicans_SC5314" = c(980040,983792),
	# 					"Ca21chr7_C_albicans_SC5314" = c(425812,428712),
	# 					"Ca21chrR_C_albicans_SC5314" = c(1743190,1747664) )

	color_1 = "#E7872B" # Orange
	color_2 = "#825CA6" # Purple
	color_3 = "#5AAA46" # Greene

	univ_pos = absolute_positions( chr_lengths )

	suppressWarnings({
		if ( graphs == "g" ){
			geno_graphs( path_and_file, chr_labs, chr_lengths, univ_pos )
			recom_graphs( path_and_file, chr_labs, chr_lengths, univ_pos )
		}

		if ( graphs == "s" ){
			shared_geno_graphs( path_and_file, chr_labs, chr_lengths, univ_pos )
		}

		if ( graphs == "c" ) {
			centro_ymax = 600 # This will depend on the number of progeny you have; remember this includes recombinations across all of them.  
			model_dat = modelling_data_prep( path_and_file, chr_lengths )
			centro_graphs( path_and_file, centromeres, chr_lengths, univ_pos, chr_labs, model_dat, ymax = centro_ymax)
		}

		if ( graphs == "h" ) {
			rhg_filename = args[4]

			recom_histo( rhg_filename, path_and_file, bandwidth = 5 )
			track_lengths( path_and_file, chr_labs, chr_lengths, univ_pos, bandwidth = 10 )
		}
	})

}


main()
warnings()