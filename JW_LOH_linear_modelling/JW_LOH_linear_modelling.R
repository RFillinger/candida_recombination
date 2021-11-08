library("ggplot2")

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

make_the_LOH_graph <- function( path_and_file, chr_labs, chr_lengths, centromeres ){

	df <- read.csv( path_and_file )

	histo <- ggplot(df, aes( x = df[,3], y = `recoms` )) + 
		geom_bar( stat = "identity" ) +
		lims( y = c(0,100) ) +
		labs( y = "Frequency", x = "Recombination events") 

	ggsave( paste0(path_and_file,"_recom_histo.png"), width = 7, height = 7)
	
	unique_chr = unique( df[,1])

	arm <- c()
	
	for (chromosomes in unique_chr ){

		chr_subset <- df[ df$"chr" == chromosomes, ]
		
		for ( row in 1:nrow(chr_subset) ){

			if ( chr_subset[row,2] < centromeres[[chromosomes]][1] ){
				arm <- c(arm, "L")

			} else if ( chr_subset[row,2] > centromeres[[chromosomes]][2] ) {
				arm <- c(arm, "R")

			} else {
				arm <- c(arm, "C")

			}

		}
	
	}

	df <- cbind( df, arm )

	for (chromosomes in unique_chr ){

		subset = df[ df[,1] == chromosomes, ]
		left_arm = subset[ subset$"arm" == "L",]
		right_arm = subset[ subset$"arm" == "R",] 

		left_formula = left_arm$"recoms" ~ left_arm[,3]
		right_formula = right_arm$"recoms" ~ right_arm[,3]

		left_mod = lm( left_formula )
		right_mod = lm( right_formula )

		right_r2_val = round(summary(right_mod)$r.squared, 2 )
		left_r2_val = round(summary(left_mod)$r.squared, 2 )

		right_pval = lmp( right_mod )
		left_pval = lmp( left_mod )

		print( paste0( "Chromosome: ", chromosomes ))
		print( paste0( "Right arm: P-val: ", right_pval, " R-squared: ", right_r2_val ))
		print( paste0( "Left arm: P-val: ", left_pval, " R-squared: ", left_r2_val ))

	}

}

make_the_LOH_graph_2 <- function( path_and_file, chr_labs, chr_lengths, centromeres ){

	df <- read.csv( path_and_file )

	df <- df[,c(1:28)]

	# Make a new column in the df to show the 100-%Heterozygous ( the number of homozygous e.g. the number of LOH events )

	df$"hom_tally" <- 21 - rowSums(df[,c(6:26)])/20
	
	histo <- ggplot(df, aes( x = df[,2], y = `hom_tally` )) + 
		geom_bar( stat = "identity" ) +
		lims( y = c(0,25) ) +
		labs( y = "Frequency", x = "Recombination events") 

	ggsave( "histo.png", width = 7, height = 7)
	
	unique_chr = unique( df[,3])

	arm <- c()
	
	for (chromosomes in unique_chr ){

		chr_subset <- df[ df$"Chr.." == chromosomes, ]
		
		for ( row in 1:nrow(chr_subset) ){

			left_arm_spots = c()
			left_arm_vals = c()

			right_arm_spots = c()
			right_arm_vals = c()

			names <- c( "pos", "homs" )

			ra_df = setNames(data.frame(matrix(ncol = 2, nrow = 0)), names)
			la_df = setNames(data.frame(matrix(ncol = 2, nrow = 0)), names)

			centro = c()

			if (chromosomes == 8){

				chr_name = "Ca21chrR_C_albicans_SC5314"

			} else{
			 
				chr_name = paste0( "Ca21chr", chromosomes, "_C_albicans_SC5314" )

			}	

			if ( chr_subset[row,4] < centromeres[[chr_name]][1] ){
				# print( c(chr_subset[row,2], chr_subset[row,29]) )
				left_arm_spots = c(left_arm_spots, chr_subset[row,2] )
				left_arm_vals = c( left_arm_vals, chr_subset[row,29] )
				la_df = la_df[nrow(la_df)+1, ] = c(chr_subset[row,2], chr_subset[row,29])
				arm <- c(arm, "L")

			} else if ( chr_subset[row,4] > centromeres[[chr_name]][2] ) {

				right_arm_spots = c( right_arm_spots, chr_subset[row,2] )
				right_arm_vals = c( right_arm_vals, chr_subset[row,29] )
				ra_df = ra_df[nrow(ra_df)+1, ] = c(chr_subset[row,2], chr_subset[row,29])
				arm <- c(arm, "R")

			} else {

				centro = c(centro, c(chr_subset[row,2], chr_subset[row,29]) )
				arm <- c(arm, "C")

			}

		}
	
	}

	df <- cbind( df, arm )

	for (chromosomes in unique_chr ){

		subset = df[ df[,3] == chromosomes, ]
		left_arm = subset[ subset[,30] == "L",]
		right_arm = subset[ subset[,30] == "R",] 

		left_formula = left_arm$"hom_tally" ~ left_arm[,2]
		right_formula = right_arm$"hom_tally" ~ right_arm[,2]

		left_mod = lm( left_formula )
		right_mod = lm( right_formula )

		right_r2_val = round(summary(right_mod)$r.squared, 2 )
		left_r2_val = round(summary(left_mod)$r.squared, 2 )

		right_pval = lmp( right_mod )
		left_pval = lmp( left_mod )

		print( paste0( "Chromosome: ", chromosomes ))
		print( paste0( "Right arm: P-val: ", right_pval, " R-squared: ", right_r2_val ))
		print( paste0( "Left arm: P-val: ", left_pval, " R-squared: ", left_r2_val ))

	}

}

main <- function() {

	# args = commandArgs( trailingOnly = TRUE )
	# path_and_file = args[1]

	
	path_and_file_1 = "recom_529L_SC5314_4way.csv"
	path_and_file_2 = "recom_P60002xSC_four_way.csv"
	path_and_file_3 = "Figure_4_Het_Density_Heat_Map_Data.csv"

	chr_labs = list("Ca21chr1_C_albicans_SC5314" = "1", 
					"Ca21chr2_C_albicans_SC5314" = "2",
					"Ca21chr3_C_albicans_SC5314" = "3",
					"Ca21chr4_C_albicans_SC5314" = "4",
					"Ca21chr5_C_albicans_SC5314" = "5",
					"Ca21chr6_C_albicans_SC5314" = "6",
					"Ca21chr7_C_albicans_SC5314" = "7",
					"Ca21chrR_C_albicans_SC5314" = "R" )

	chr_lengths = list( "Ca21chr1_C_albicans_SC5314" = 3190000,
						"Ca21chr2_C_albicans_SC5314" = 2230000,
						"Ca21chr3_C_albicans_SC5314" = 1800000,
						"Ca21chr4_C_albicans_SC5314" = 1600000,
						"Ca21chr5_C_albicans_SC5314" = 1190000,
						"Ca21chr6_C_albicans_SC5314" = 1030000,
						"Ca21chr7_C_albicans_SC5314" = 950000,
						"Ca21chrR_C_albicans_SC5314" = 2290000 )

	centromeres = list( "Ca21chr1_C_albicans_SC5314" = c(1563038,1565967), 
						"Ca21chr2_C_albicans_SC5314" = c(1927255,1930214),
						"Ca21chr3_C_albicans_SC5314" = c(823333,826481),
						"Ca21chr4_C_albicans_SC5314" = c(992579,996216),
						"Ca21chr5_C_albicans_SC5314" = c(468716,471745),
						"Ca21chr6_C_albicans_SC5314" = c(980040,983792),
						"Ca21chr7_C_albicans_SC5314" = c(425812,428712),
						"Ca21chrR_C_albicans_SC5314" = c(1743190,1747664) )
	print( " " )
	print( path_and_file_1 )
	make_the_LOH_graph( path_and_file_1, chr_labs, chr_lengths, centromeres )
	print( " " )
	print( path_and_file_2 )
	make_the_LOH_graph( path_and_file_2, chr_labs, chr_lengths, centromeres )
	print( " " ) 
	print( path_and_file_3 )
	make_the_LOH_graph_2( path_and_file_3, chr_labs, chr_lengths, centromeres )


	}



main()
