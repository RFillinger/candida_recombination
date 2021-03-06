
library("stringr")
library("ggplot2")
library("grid")
library("chisq.posthoc.test")
library("ggpubr")
library("diptest")

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
							ymax = 1, tick_perct = 0.025, add_stdev_lines = 0, color = "gold", theme_color = "black" ) {

	path = path_components( path_and_file )[1]
	file = path_components( path_and_file )[2]

	df = read.csv( path_and_file )

	# Chi-squared test (goodness of fit) lines go here
	# equal_probability_vector = rep((1/length(df$recom.by.mkr)), length(df$recom.by.mkr)) # Assumes random chance of a recombination for each bin. 
	# redistribution_vector = rep( round(sum(df$recom.by.mkr)/length(df$recom.by.mkr), 0), length(df$recom.by.mkr) )

	# chi_squared = chisq.test( x = df$recom.by.mkr, p = equal_probability_vector )
	# chi_squared = chisq.test( x = df$recom.by.mkr, y= redistribution_vector )
	# print( chi_squared$expected  ) # Every value should be above 5
	# print( chi_squared )

	# shap_test = shapiro.test(df$recom.by.mkr) # Test for normalcy of distributions
	# levene.test( df$recom.by.mkr )

	recom_stdev = sd(df$recom.by.mkr)
	recom_avg = mean(df$recom.by.mkr)
	threshold = recom_avg + (2*recom_stdev) # represents 95% of data

	big_bins = df[df$recom.by.mkr>threshold,]

	y_pos = -0.5
	col_1 = "magenta"
	col_2 = "turquoise"

	# histogram <- ggplot(df, aes(x=df$recom.by.mkr)) + 
	# 	geom_histogram(binwidth=1) + 
	# 	lims( x=c(0,150) ) +
	# 	labs( y = "Frequency", x = "Recombination events") + 
	# 	geom_segment( x = recom_avg, y = 0, xend = recom_avg, yend = y_pos, color = col_1) + 
	# 	geom_segment( x = threshold, y = 0, xend = threshold, yend = y_pos, color = col_2 ) + 
	# 	geom_segment( x = opp_threshold, y = 0, xend = opp_threshold, yend = y_pos, color = col_2 ) + 
	# 	annotation_custom( textGrob("Mean", gp = gpar( fontsize = 10, col = col_1)), # Adds label to mean line
	# 						xmin = recom_avg , xmax = recom_avg, ymin = y_pos, ymax = y_pos ) + 
	# 	annotation_custom( textGrob("Mean + 2σ", gp = gpar( fontsize = 10, col = col_2 )), # Adds label to mean line
	# 						xmin = threshold , xmax = threshold, ymin = y_pos, ymax = y_pos )

	# ggsave( paste0(path, substr( file, 1, str_length(file)-4 ), ".recom_hist", file_extension ), width=7, height=7 )

	chr_vec = c()
	chr_reco_vec = c()

	chromosomes <- unique( df$chr )
	max_sizes = vector( mode="list", length=length( chromosomes ) )
	names( max_sizes ) = as.vector(chromosomes)
	for ( chromosome in chromosomes ){
		# Subset the data frame to obtain the maximum sizes of chromosomes in the df
		chr_subset = df[df$chr == chromosome, ]
		max_sizes[[chromosome]] <- max(chr_subset$univ.pos) 

		# Find the number of rows (bins) in the chromosome and the sum of all normalized recombinations
		bin_num = nrow( chr_subset )
		sum_totes = sum(chr_subset$recom.by.mkr)
		chromosomal_recom = sum_totes/bin_num
		
		chr_vec = c(chr_vec, chromosome)
		chr_reco_vec = c(chr_reco_vec, chromosomal_recom)

	}

	cnr = data.frame( chr_vec, chr_reco_vec )
	v <- ggplot( data = cnr, aes( x = chr_vec, y = chr_reco_vec ) ) + 
		geom_col( color = "cyan" ) + 
		theme(  axis.title.y = element_blank(),
				axis.title.x = element_blank(),
				axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
		lims(y = c(0,0.06))
	ggsave( paste0(path, substr( file, 1, str_length(file)-4 ), "_NORMALIZED_RECOMS", file_extension), width=7, height=7 )


	chromosomes_delineations = unlist(max_sizes)[c(1:length(chromosomes)-1)] 
	y_disp = -(ymax*tick_perct) + 15

	p <- ggplot( data = df, aes( x = `univ.pos`, y = `recom.by.mkr` )  ) + 
		geom_col( color = "cyan" ) + 
		# geom_vline( xintercept = chromosomes_delineations, size = 0.2, color = color ) + 
		labs( y = "Recombinations per Marker", x = NULL ) + 
		ylim( c(0,ymax)) + 
		theme_bw() + 
		theme( 
		 	panel.grid.major.x = element_blank(), 
			panel.grid.minor.x = element_blank(), 
			panel.grid.minor.y = element_blank(),
			panel.border = element_blank(),
			# axis.title.x = element_blank(), 
			axis.ticks.x = element_blank(), 
			axis.text.x = element_blank(),
			# axis.title.y = element_blank(), 
			# axis.text.y = element_blank(), 
			# axis.ticks.y = element_blank()
			)

	if ( add_stdev_lines ){
		p <- p + geom_hline( yintercept = recom_avg ) + 
				 geom_hline( yintercept = threshold )

		p <- p + 
		geom_segment( x = min(unlist(max_sizes)), y = recom_avg, xend = max(unlist(max_sizes)), yend = recom_avg, color = col_1) + 
		geom_segment( x = min(unlist(max_sizes)), y = threshold, xend = max(unlist(max_sizes)), yend = threshold, color = col_2) + 
		annotation_custom( textGrob("Mean", gp = gpar( fontsize = 10, col = col_1)), # Adds label to mean line
							xmin = 0 , xmax = 0, ymin = recom_avg + (ymax*tick_perct), ymax = recom_avg + (ymax*tick_perct) ) + 
		annotation_custom( textGrob("Mean + 2σ", gp = gpar( fontsize = 10, col = col_2 )), # Adds label to mean line
							xmin = 0 , xmax = 0, ymin = threshold + (ymax*tick_perct), ymax = threshold + (ymax*tick_perct) )
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
		p <- p + annotation_custom( textGrob(chr_labs[[chromosomes]], gp = gpar( fontsize = 17, col = theme_color)), # Adds the chromosome labels
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
							ymax = 200, tick_perct = 0.01, theme_color = "black", line_Width = .1 ) {

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

	xmax = max(df$absolute_pos)

	p <- ggplot( data = df, aes( x = `absolute_pos`, y = `cumulative_recom_cnt` )  ) + 
		geom_step( aes(color = `chromosome`, lwd = line_Width ) ) + 
		theme_bw() +  
		theme( legend.position = "none", axis.text.x = element_blank() ) + 
		ylim( c(0,ymax)) + 
		# labs( x = NULL, y = "Cumulative Recombination" )

		theme( 
		 	panel.grid.major.x = element_blank(), 
		 	# panel.grid.major.y = element_blank(),
			panel.grid.minor.x = element_blank(), 
			panel.grid.minor.y = element_blank(),
			panel.border = element_blank(),
			axis.title.x = element_blank(), 
			axis.ticks.x = element_blank(), 
			axis.text.x  = element_blank(),
			axis.title.y = element_blank(), 
			axis.text.y  = element_blank(), 
			axis.ticks.y = element_blank()
			) + 

		geom_segment( x = 0, y = 0,      xend = 0,    yend = ymax,   color = theme_color) + # y-Axis
		geom_segment( x = 0, y = 0,      xend = xmax, yend = 0,      color = theme_color) + # x-Axis
		geom_segment( x = 0, y = ymax/2, xend = xmax, yend = ymax/2, color = theme_color) + # Half way point
		geom_segment( x = 0, y = ymax,   xend = xmax, yend = ymax,   color = theme_color)   # Tippy top
		
	# I hate that I have to do this... it adds ticks and lines and editable labels to the graph where the centromeres are. Ugh.
	y_disp = -(ymax*tick_perct) - 5 

	i = 1 
	for ( centromere in centro_means ){
		p <- p + geom_segment( x = centromere, y = 0, xend = centromere, yend = tick_perct * ymax, color = theme_color ) # Adds the better custom ticks
		# p <- p + geom_segment( x = centromere, y = 0, xend = centromere, yend = tick_perct * ymax, color = theme_color ) +
		# 		annotation_custom( textGrob(chr_labs[[i]], gp = gpar( fontsize = 17, col = theme_color)), # Adds the chromosome labels
		# 			xmin = centromere , xmax = centromere, ymin = y_disp, ymax = y_disp ) 
		i = i + 1
	}

	ggsave( paste0(path, substr( file, 1, str_length(file)-4 ), file_extension), width=7, height=7 )

	# i = 0
	# for (arms in names(model_dat)) {

	# 	chr_name = unique( model_dat[[ arms ]]$chromosome )

	# 	if (i%%2 == 1){
	# 		x_loc = round(chr_lengths[[chr_name]]*(1/4)) + universal_positions[[chr_name]]
	# 		y_loc = ymax*0.9

	# 	} else {
	# 		x_loc = round(chr_lengths[[chr_name]]*(3/4)) + universal_positions[[chr_name]]
	# 		y_loc = ymax*0.7
	# 	}

	# 	arm_model = lm (model_dat[[ arms ]]$cumulative_recom_cnt ~ model_dat[[ arms ]]$absolute_pos)
	# 	r2_val = round(summary(arm_model)$r.squared, 2 )
	# 	model_p_val = lmp( arm_model )
	# 	coeff_p_val = round( summary( arm_model)$coefficients[,4], 2)

	# 	lm_model = model_dat[[arms]][,4] ~ model_dat[[arms]][,3]
		
	# 	# Prints p-values of model fit to each arm
	# 	# print( arms )
	# 	# print( model_p_val )
		
	# 	# This set will generate lines, but make a lot of warnings
	# 	# p = p + geom_smooth( data = model_dat[[arms]], method = "lm", se = FALSE, color = "black", na.rm = TRUE ) + 
	# 	# 		annotation_custom( textGrob( r2_val, gp = gpar( fontsize = 7, col = theme_color)), # Adds the r squared values
	# 	# 			xmin = x_loc , xmax = x_loc, ymin = y_loc, ymax = y_loc )

	# 	# This chunk won't make lines or warnings
	# 	p = p + geom_smooth( data = model_dat[[arms]], formula = lm_model, method = "lm", se = FALSE, color = "black", na.rm = TRUE ) + 
	# 			annotation_custom( textGrob( r2_val, gp = gpar( fontsize = 7, col = theme_color)), # Adds the r squared values
	# 				xmin = x_loc , xmax = x_loc, ymin = y_loc, ymax = y_loc ) #+ 
	# 			# annotation_custom( textGrob( model_p_val, gp = gpar( fontsize = 7, col = theme_color)), # Adds the p-values
	# 			# 	xmin = x_loc, xmax = x_loc, ymin = y_loc - 0.05*ymax, ymax = y_loc - 0.05*ymax )
	# 	i = i + 1
	# }

	# ggsave( paste0(path, substr( file, 1, str_length(file)-4 ), "_lm", file_extension ), width=7, height=7 )

}


shared_geno_graphs <- function( path_and_file, chr_labs, chr_lengths, univ_pos, a = "1", b = "2", h = "n", missing = "-", include_missing = 0, 
								color_1 = "#E7872B", color_2 = "#825CA6", color_3 = "#5AAA46" ){
	# color_1 = "#E7872B" # Orange  color_2 = "#825CA6" # Purple # color_3 = "#5AAA46" # Greene
	path = path_components( path_and_file )[1]
	file = path_components( path_and_file )[2]
	df = read.csv( path_and_file )

	if (include_missing == 0 ){

		df = df[ df$allele != missing, ]
		
		df$allele = factor( df$allele, levels = c(b,h,a) )
		colors = c( color_2, color_3, color_1 )

		df = df[ order(df$a_prop), ]
		df$Progeny.Strain = factor( df$Progeny.Strain, levels = unique(df$Progeny.Strain) )
		
	}

	p <- ggplot( data = df, aes( fill = `allele`, x = `Progeny.Strain`, y = `prop` )  ) +
			geom_bar( position = "fill", stat = "identity", width = 2.5 ) + 
			scale_fill_manual( values = colors ) + 
			scale_x_discrete( limits = df$Progeny.Strain) + # Adds parent and progeny names 
			theme_bw() + 
			theme(  axis.text = element_blank(), 
					axis.title = element_blank(), 
					legend.text = element_blank(),
					axis.text.x = element_blank(),
					panel.border = element_blank(),
					panel.grid.major = element_blank(),
					panel.grid.minor = element_blank(), 
					axis.ticks = element_blank(), 
					legend.position = "none" )
			
	ggsave( paste0(path, substr( file, 1, str_length(file)-4 ), file_extension ), width=10, height=7)

}


track_lengths <- function( path_and_file, chr_labs, chr_lengths, univ_pos, bwidth = 0.1, ymaximum = 500, xmaximum = 750, theme_color = "black" ){

	path = path_components( path_and_file )[1]
	file = path_components( path_and_file )[2]
	df = read.csv( path_and_file )

	# print(head(df))

	g <- ggplot( data = df, aes( x = `Recombination.Event.Length..bp.`/1000)) + # kilo base pairs
	geom_histogram( binwidth = bwidth, color = "black", fill = "white") +
	labs( x = paste0("Estimated Recombination Size (Kbp) Binsize: ", bwidth), y = "Frequency") +
	ylim( c(0,ymaximum)) +
	xlim( c(0,xmaximum)) + 

	theme_bw() + 
	theme( 
		 	panel.grid.major.x = element_blank(), 
			panel.grid.minor.x = element_blank(), 
			panel.grid.minor.y = element_blank(),
			panel.border = element_blank(),
			# axis.title.x = element_blank(), 
			axis.ticks.x = element_blank(), 
			axis.text.x = element_blank(),
			# axis.title.y = element_blank(), 
			# axis.text.y = element_blank(), 
			# axis.ticks.y = element_blank()
			)

	# Loop for the delineations
	# x_spots = c(0,500,1000,1500,2000)
	x_spots = c(0,250,500,750)

	for ( tickmarks in x_spots ){
		g <- g + 
			geom_segment( x = tickmarks, y = 0, xend = tickmarks, yend = -0.015* ymaximum, color = theme_color) +
			annotation_custom( textGrob(as.character(tickmarks), gp = gpar( fontsize = 17, col = theme_color )), # Adds the chromosome labels
						xmin = tickmarks , xmax = tickmarks, ymin = -0.035* ymaximum, ymax = -0.035* ymaximum )
	}

	ggsave( paste0( path, substr( file, 1, str_length(file)-4 ), file_extension ), width=7, height=7 )

}


recom_histo <- function( rhf_filename, path_and_file, bwidth, ymaximum = 35, theme_color = "black" ) {

	path_stuff = path_components( path_and_file )
	path = path_stuff[1]
	file = path_stuff[2]

	df = read.csv( rhf_filename )

	recombinations = df[,2]

	x <- ggplot( data = df, aes( x = recombinations )) + 
	geom_histogram( binwidth = bwidth, color = "black", fill = "white" ) + 
	labs( x = "", y = "Recombination Frequency" ) + 
	theme_bw() + 
	ylim( c(0,ymaximum)) +
	xlim(c(0,105)) + 
	# xlim( c(0,100)) + 
	theme( 
		 	panel.grid.major.x = element_blank(), 
			panel.grid.minor.x = element_blank(), 
			panel.grid.minor.y = element_blank(),
			panel.border = element_blank(),
			# axis.title.x = element_blank(), 
			axis.ticks.x = element_blank(), 
			axis.text.x = element_blank(),
			# axis.title.y = element_blank(), 
			# axis.text.y = element_blank(), 
			# axis.ticks.y = element_blank()
			)

	# Loop for the delineations
	x_spots = c(5,25,50,75,100)
	# x_spots = c(0,100,200,300,400,500,600)
	for ( tickmarks in x_spots ){
		x <- x + geom_segment( x = tickmarks, y = 0, xend = tickmarks, yend = -0.015* ymaximum, color = theme_color) +
				annotation_custom( textGrob(as.character(tickmarks), gp = gpar( fontsize = 17, col = theme_color )), # Adds the chromosome labels
						xmin = tickmarks , xmax = tickmarks, ymin = -0.035* ymaximum, ymax = -0.035* ymaximum )
	}

	ggsave( paste0( path, "RECO_histo_", substr(file,12,str_length(file)-4 ), ".png" ), width=7, height=7 )

	# dt = dip.test(recombinations)
	# print(paste0("Results of dip test for modality:"))
	# print(paste0("Test statistic: ", dt$statistic))
	# print(paste0("P value: ", dt$p.value))

}


loh_reco <- function( path_and_file, chr_labs, chr_lengths, universal_positions,  
							ymax = 1, tick_perct = 0.025, theme_color = "black" ){ 

	path = path_components( path_and_file )[1]
	file = path_components( path_and_file )[2]

	df = read.csv( path_and_file )

	recom_stdev = sd(df$recom.by.mkr)
	recom_avg = mean(df$recom.by.mkr)
	threshold = recom_avg + (2*recom_stdev) # represents 95% of data

	big_bins = df[df$recom.by.mkr>threshold,]

	y_pos = -0.5

	chr_vec = c()
	chr_reco_vec = c()

	chromosomes <- unique( df$chr )
	max_sizes = vector( mode="list", length=length( chromosomes ) )
	names( max_sizes ) = as.vector(chromosomes)
	for ( chromosome in chromosomes ){
		# Subset the data frame to obtain the maximum sizes of chromosomes in the df
		chr_subset = df[df$chr == chromosome, ]
		max_sizes[[chromosome]] <- max(chr_subset$univ.pos) 

		# Find the number of rows (bins) in the chromosome and the sum of all normalized recombinations
		bin_num = nrow( chr_subset )
		sum_totes = sum(chr_subset$recom.by.mkr)
		chromosomal_recom = sum_totes/bin_num
		
		chr_vec = c(chr_vec, chromosome)
		chr_reco_vec = c(chr_reco_vec, chromosomal_recom)

	}


	chromosomes_delineations = unlist(max_sizes)[c(1:length(chromosomes)-1)] 
	y_disp = -(ymax*tick_perct) + 15

	p <- ggplot( data = df, aes( x = `univ.pos`, y = `recom.by.mkr` )  ) + 
		geom_col( ) + 
		# geom_vline( xintercept = chromosomes_delineations, size = 0.2, color = color ) + 
		labs( y = "Recombinations per Marker", x = NULL ) + 
		ylim( c(0,ymax)) + 
		theme_bw() + 
		theme( 
		 	panel.grid.major.x = element_blank(), 
			panel.grid.minor.x = element_blank(), 
			panel.grid.minor.y = element_blank(),
			panel.border = element_blank(),
			axis.title.x = element_blank(), 
			axis.ticks.x = element_blank(), 
			axis.text.x = element_blank(),
			axis.title.y = element_blank(), 
			axis.text.y = element_blank(), 
			axis.ticks.y = element_blank()
			)

	# One loop for the delineations 
	for ( chromosomes in names(chromosomes_delineations) ){
		x_spot = chr_lengths[[chromosomes]] + universal_positions[[ chromosomes ]]
		p <- p + geom_segment( x = x_spot, y = 0, xend = x_spot, yend = -tick_perct * ymax, color = theme_color)
	}
	

	# Second loop for the labels that go under the delineations
	i = 1 
	for ( chromosomes in names(chr_labs) ){
		x_spot = round(chr_lengths[[ chromosomes ]]/2 + universal_positions[[ chromosomes ]])
		p <- p + annotation_custom( textGrob(chr_labs[[chromosomes]], gp = gpar( fontsize = 17, col = theme_color)), # Adds the chromosome labels
						xmin = x_spot , xmax = x_spot, ymin = -tick_perct * ymax, ymax = -tick_perct * ymax )
		i = i + 1
	}

	ggsave( paste0(path, "LOH_", substr( file, 1, str_length(file)-4 ), file_extension), width=7, height=7 )

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

	color_1 = "#E7872B" # Orange
	color_2 = "#825CA6" # Purple
	color_3 = "#5AAA46" # Greene

	color_4 = "#EBEDEE" # Gray from the theme_bw()

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
			centro_ymax = 500 # This will depend on the number of progeny you have; remember this includes recombinations across all of them.  
			model_dat = modelling_data_prep( path_and_file, chr_lengths )
			centro_graphs( path_and_file, centromeres, chr_lengths, univ_pos, chr_labs, model_dat, ymax = centro_ymax)
		}

		if ( graphs == "h" ) {
			rhg_filename = args[4]

			recom_histo( rhg_filename, path_and_file, bwidth = 5 )
			track_lengths( path_and_file, chr_labs, chr_lengths, univ_pos, bwidth = 10 )
		}

		if ( graphs == "loh" ){ 

			loh_reco( path_and_file, chr_labs, chr_lengths, univ_pos )

		}
	})

}


main()
warnings()