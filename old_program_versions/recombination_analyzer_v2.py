
import statistics
import operator
import random
import numpy
import copy
import sys
import os

import time
import pdb

def csv_reader( input_file ):

	'''Reads in .csv files into a list'''

	line_list = []
	for lines in input_file: 
		line_list.append( lines.strip().split(",") )

	return line_list 

def transpose(l1):

	""" Transposes a 2D python array """

	new_list = []
	# iterate over list l1 to the length of an item
	for i in range(len(l1[0])):
		# print(i)
		row =[]
		for item in l1:
			# appending to new list with values and index positions
			# i contains index position and item contains values
			row.append(item[i])
		new_list.append(row)

	return new_list


def csv_printer( csv_list, output_file ):

	for items in csv_list: 
		line = ",".join(items)
		print( line, file = output_file)


def parse_path( path_n_name, sep = "/" ): 

	''' Takes a file path and name and splits it into two variables to be used for output. '''

	path_list = path_n_name.split( sep )
	file_name = path_list[-1]
	dir_path = "/".join(path_list[:-1]) + "/"
	if dir_path == "/": 
		dir_path = ""

	return dir_path, file_name



def albicans_header_creator(length): 
	
	chr_nums = [ 1, 2, 3, 4, 5, 6, 7, "R" ]
	chr_lengths = [ 3190000, 2230000 , 1800000, 1600000, 1190000, 1030000, 950000, 2290000 ]

	header = []
	chr_ID = [ ]
	catalog_ID = []
	chr_loc = []
	
	for i in range( 1, length+1 ):
		catalog_ID.append( i ) # Just enumerate this line, so it can be used. 
		num = random.randint(1,8)
		chr_ID.append( num )

	chr_ID.sort()

	for n, elements in enumerate( chr_ID ): 

		if elements == 8: 

			chr_ID[n] = "R"

	for index, chromosome in enumerate(chr_nums):

		cnt = chr_ID.count( chr_nums[index] )
		locations = random.sample( range(1, chr_lengths[index]), cnt )
		
		try:

			locations.sort()

			for things in locations: 
				chr_loc.append( things )

		except TypeError: 
			print( "ERROR" )
			print( "cnt: ", cnt )
			print( "index: ", index )
			print( chr_nums )

	for i, chromosomes in enumerate(chr_ID):

		chr_ID[i] = "Ca21chr" + str(chr_ID[i]) + "_C_albicans_SC5314"
		
	
	chr_ID.insert( 0, "Chromosome" )
	chr_loc.insert( 0, "Chr_Position" )
	catalog_ID.insert( 0, "# Catalog ID" )

	header = [ chr_ID, chr_loc, catalog_ID ]
	
	return header


def test_creator( output_file_name, number_of_prog_to_generate = 66, length = 200, inverse_missing = 4 ): 
	""" Creates a 4-way allele file for testing the other components of this program. """

	shoot_trouble = 0
	
	parent_combs = { 1: ["aa", "aa"], 
					 2: ["aa", "ab"], 
					 3: ["aa", "bb"], 
					 4: ["ab", "aa"],
					 5: ["ab", "ab"],
					 6: ["ab", "bb"],
					 7: ["aa", "bc"],
					 8: ["ab", "ac"],
					 9: ["ab", "bc"],
					10: ["ab", "cc"],
					11: ["ab", "cd"] }

	prog_combs = ["aa", "ab", "bb", "ac", "bc", "cc", "ad", "bd", "cd", "dd"]

	# Make a file containing all combinations of alleles in a line and see where the recombinations are. 
	if shoot_trouble: 
		dist_list = [0] * 11 # creates a list showing the distribution of random numbers

	big_boi = []

	header = [ "Parent 1", "Parent 2" ]
	for x in range( 0, number_of_prog_to_generate ): 
		string = "Progeny " + str( x )
		header.append( string )

	big_boi.append( header )

	index = 1

	while index <= length: 

		num = random.randint( 1, 11 )

		if shoot_trouble: 
			dist_list[ num - 1 ] += 1

		leel_list = []

		leel_list = copy.deepcopy( parent_combs[num] ) # Deepcopy because I'm worried about linking the list to the dictionary

		for x in range( 0, number_of_prog_to_generate ):

			missing = 0
			missing_key = random.randint( 1, inverse_missing ) #Introduces missing data points as "-" just like in real data. 
			if missing_key == 1: 
				missing = 1

			if (num == 1) and not missing: # only "aa" is available
				
				leel_list.append( prog_combs[0] ) 

			elif ( num > 1 ) and ( num <= 6 ) and not missing: # Only "a" or "b" allele combinations are possible

				rint_6 = random.randint(1, 3)
				leel_list.append( prog_combs[rint_6-1])

			elif ( num > 6 ) and ( num <= 10 ) and not missing: # A, B, and C alleles are available

				rint_10 = random.randint(1, 6)
				leel_list.append( prog_combs[rint_10-1])

			elif (num == 11) and not missing: # All aleles A through D are available

				rint_11 = random.randint(1, 10)
				leel_list.append( prog_combs[rint_11-1])

			else:
				leel_list.append( "-")

		big_boi.append( leel_list )

		index += 1

	if shoot_trouble: 
		print( big_boi )
		print( dist_list )
		print( type( big_boi) )

	trp_big_boi = transpose( big_boi )

	header_parts = albicans_header_creator(length)

	type_swap_header_list = [ [], [], [] ]
	
	for n, sublists in enumerate( header_parts ): 

		for elements in sublists: 

			type_swap_header_list[n].append( str(elements) )


	test_file = open( "test/test.csv", "w" )

	print( ",".join(type_swap_header_list[0]), file = test_file )
	print( ",".join(type_swap_header_list[1]), file = test_file )
	print( ",".join(type_swap_header_list[2]), file = test_file )


	csv_printer( trp_big_boi, test_file ) # Change this to print the file and we're good to go. 
	test_file.close()


def marker_cleaner(four_way_file_name): 
	"""This function removes markers that are missing from eiter parent and markers that are
	ambiguous between parents (the parents share an allele). """

	print_undesirables = 1

	dir_path, file_name = parse_path( four_way_file_name )

	four_way_file = open( four_way_file_name, "r" )

	untrp_leels = csv_reader( four_way_file ) # Load the file, but it needs to be transposed
	leels = transpose( untrp_leels )

	keepers = []
	chuckers = []

	header = 1

	for lines in leels: 

		if header: 
			keepers.append( lines )
			chuckers.append( lines )
			header = 0
			continue

		parent1 = lines[3]
		parent2 = lines[4]
		
		try:
			lines[1] = int( lines[1] ) # Convert to ints so they can get sorted properly
			lines[2] = int( lines[2] )
		except IndexError: 
			print( "Something is wrong with the format in this entry: ", lines )
			continue

		if ("-" in parent1) or ("-" in parent2): # Remove markers missing from parents
			chuckers.append( lines )
			continue

		if ("a" in parent1) and ("a" in parent2): # Remove ambiguous markers; only need to check for a and b 
			chuckers.append( lines )
			continue

		if ("b" in parent1) and ("b" in parent2): # Remove ambiguous markers 
			chuckers.append( lines )
			continue

		keepers.append( lines )

	print( "marker_cleaner removed ", len(chuckers) - 1 , " markers (", 100*(len(chuckers)-1)/(len(chuckers)-1+len(keepers)-1), "%","removed )" )

	srt_keepers = sorted( keepers, key=operator.itemgetter(0,1) )
	for lines in srt_keepers: 
		lines[1] = str(lines[1]) # Convert back to strings so they can be output. 
		lines[2] = str(lines[2])
	
	header = srt_keepers.pop(-1) # Move the header string back to the front
	srt_keepers.insert(0, header)

	trp_keepers = transpose( srt_keepers ) # File needs to be tranposed again. 
	new_file = open( dir_path + "cleaned_" + file_name, "w" )
	for lines in trp_keepers: 
		print( ",".join(lines), file = new_file )
	new_file.close()

	if print_undesirables:

		srt_chuckers = sorted( chuckers, key=operator.itemgetter(0,1) )
		for lines in srt_chuckers: 
			lines[1] = str(lines[1]) # Convert back to strings so they can be output. 
			lines[2] = str(lines[2]) 
		
		header = srt_chuckers.pop(-1) # Move the header string back to the front
		srt_chuckers.insert(0, header)

		trp_chuckers = transpose( srt_chuckers )
		new_un_file = open( dir_path + "undesirables_" + file_name, "w" )
		for lines in trp_chuckers: 
			print( ",".join(lines), file = new_un_file )
		new_un_file.close()


def four_to_F2( cleaned_file_name ): 
	
	"""Converts 4-way alele data to F2 
	ALL INPUT INTO THIS FUNCTION MUST GO THROUGH THE 'marker_cleaner' FUNCTION FIRST!!!"""
	dir_path, file_name = parse_path( cleaned_file_name )

	cleaned_marker_file = open( dir_path + file_name, "r" )

	untrp_cleaned_markers = csv_reader( cleaned_marker_file )

	cleaned_markers = transpose( untrp_cleaned_markers )

	f2_array = []

	header = 1
	for lines in cleaned_markers: 
		
		if header: 
			f2_array.append( lines )
			header = 0
			continue		

		parent1 = lines[3]
		parent2 = lines[4]

		new_line = copy.deepcopy(lines[:3]) # Creates the new marker for conversion and adds alleles for parent 1 and 2 (since we know they're already mutually exclusive)
		new_line.append( "1" )
		new_line.append( "2" )

		for n, elements in enumerate(lines[5:]): # This DOES NOT check for post-parasex LOH events ("dd" and perhaps "cc" are examples). Not sure how important those are, but they'll be called "n." 

			if "-" in elements: 
				new_line.append( "-" )
			elif elements == parent1: 
				new_line.append( "1" )
			elif elements == parent2: 
				new_line.append( "2" )
			else:
				new_line.append( "n" )

		f2_array.append( new_line )

	trp_f2 = transpose(f2_array)

	new_file = open( dir_path + "f2_" + file_name, "w")
	csv_printer( trp_f2, new_file )
	new_file.close()


def recombination_analysis( f2_file_name, output_file_name, a = "1", b = "2", h = "n" ):
	"""This function takes the cleaned f2 data and searches for recombination events, and stores them """

	#########
	## ISSUES: 
	## All solved, for now. 

	f2_file = open( f2_file_name, "r" )
	output_file = open( output_file_name, "w" )

	f2_array = csv_reader( f2_file )

	# Printing the stats file
	header_str = "Sample Name, Filled Markers, Parent 1 Markers, Parent 2 Markers, Heterozygous Markers, Missing Markers, Number of Recombination Events, Recombined Chromosome, Recombination Start Site (bp), Recombination Event Length (bp), Marker"
	print( header_str , file = output_file )

	chromosomes = f2_array[0]
	chr_pos     = f2_array[1]
	mrkr_cat_id = f2_array[2]
	parent1     = f2_array[3]
	parent2     = f2_array[4]
	children    = f2_array[5:]
	recomb_list = []

	for progeny in children: 

		recombs = { "prog_id": "", "a_leels": 0 , "b_leels": 0, "h_leels": 0, "missing": 0, "recombs": 0, "chro": {}} 
		recombs[ "prog_id" ] = progeny[0]

		first_marker = 1
		first_col = 1 
		for i, marker in enumerate(progeny):

			if first_col: # First column contains the progeny names and other labels; useless for now. 
				first_col = 0 
				continue

			missing_bool = marker.strip() == "-"
			# pdb.set_trace( header = "Progeny: " + str(progeny[0]) + "\nIndex (i): " + str( i ) + "\nMarker: " + marker )
			if first_marker and (not missing_bool): 
				stor_chr = chromosomes[i] # Storage variable for current chromosome
				stor_chr_pos = chr_pos[i] # When a recombination occurs, this is subtracted from the new length to obtain recombination size in bp
				stor_marker = marker      # Storage var for current marker (needs to be compared to the next marker to assess recombination)
				first_marker = 0
			
			elif not first_marker:

				if chromosomes[i] not in recombs["chro"]: # Need to create a new entry for new chromosomes in the recombinations 
					recombs["chro"][chromosomes[i]] = []

				same_chr_bool = ( chromosomes[i] == stor_chr ) # Make sure we're on the same chromosome (recombination doesn't occur across chromosomes)
				recomb_bool   = ( stor_marker != marker ) and ( not missing_bool ) and same_chr_bool # If they're different recombination has occured!
				
				# pdb.set_trace( header = "\nchromosome: " + str(chromosomes[i]) + \
					# "\nstor_chr: " + str(stor_chr) + \
					# "\nboolean_of_equality: " + str( chromosomes[i] == stor_chr ) )

				if ( not same_chr_bool ) and (not missing_bool): # If the chromosome has changed, start over.
					
					# pdb.set_trace( header = "Chromosome changed!" + \
					# 	"\nPrevious chromosome: " + str(stor_chr) + \
					# 	"\nCurrent chromosome: " + str(chromosomes[i]) + \
					# 	"\nOld marker " + str(stor_marker) + \
					# 	"\nNew marker " + str(marker))
					
					stor_chr = chromosomes[i]
					stor_chr_pos = chr_pos[i]
					stor_marker = marker

				elif recomb_bool: 

					# pdb.set_trace( header = "Recombination! " + "\nProgeny: " + str(progeny[0]) + \
					# 										"\nPrevious chromosome: " + str(stor_chr) + \
					# 										"\nCurrent chromosome: " + str(chromosomes[i]) + \
					# 										"\nOld marker " + str(stor_marker) + \
					# 										"\nNew marker " + str(marker))
					
					recombs[ "recombs" ] += 1
					recomb_size = int( chr_pos[i] ) - int( stor_chr_pos ) # Compute the estimated size of the recombination event and store it
					recombs["chro"][ chromosomes[i] ].append( [ int( stor_chr_pos ), recomb_size, marker ] ) # Store the putative recombination start location and estimated size
					
					stor_marker = marker # Change the stor_marker and stor_chr_pos because this is the new front boundary for a recombination event
					stor_chr_pos = chr_pos[i] 

			if marker == a: 
				recombs[ "a_leels" ] += 1 
			elif marker == b: 
				recombs[ "b_leels" ] += 1
			elif marker == h: 
				recombs[ "h_leels" ] += 1

		# Make the progeny data entry in the file: 
		progeny_str = recombs["prog_id"]
		filled_mkrs = recombs["a_leels"] + recombs["b_leels"] + recombs["h_leels"] 
		parent1_mkr = recombs["a_leels"]
		parent2_mkr = recombs["b_leels"]
		het_mkrs    = recombs["h_leels"]
		missing_mkr = recombs["missing"]
		recomb_cnt  = recombs["recombs"]

		print_str = progeny_str + "," + str(filled_mkrs) + "," + str(parent1_mkr) + "," + str(parent2_mkr) + "," + str(het_mkrs) + "," + str(missing_mkr) + "," + str(recomb_cnt)

		for chromos, recombinations in recombs["chro"].items():

			for events in recombinations: 

				event_start = events[0]
				event_length = events[1]
				event_leel = events[2]

				print( print_str + "," + str(chromos) + "," + str(event_start) + "," + str(event_length) + "," + event_leel, file = output_file )

	f2_file.close()
	output_file.close()


def recom_Stacker(chr_lengths, recom_stats_file_name, output_file_name, stack_type, resolution = 100000, a = "1", b = "2", h = "n" ): 
	"""Creates a dictionary of chromosomes and bin sizes that wll be used for histograms to depict recombination hotspot for each allele type.""" 

	#########
	## ISSUES: 
	## Something is wrong with the distances between bins. Either here, or in the R program.  

	total = 0
	univ_chr_dict = {}
	for chromosomes, lengths in chr_lengths.items(): 
		univ_chr_dict[ chromosomes ] = total
		total += int(lengths)

	recom_bins, a_geno_bins, b_geno_bins, n_geno_bins = {}, {}, {}, {}
	dict_list = [ recom_bins, a_geno_bins, b_geno_bins, n_geno_bins ]

	# Build a set of dictionaries to store the recombination events for each chromosome and a nested dictionary to contain the bins  
	for bins in dict_list: 
		for chromo, length in chr_lengths.items(): 
			chr_bin_dict = {}
			index = 0
			while index <= length: 
				# pdb.set_trace( header = "resolution: " + str(resolution) + "\nIndex: " + str( index ) + "\nChromosome: " + chromo )
				chr_bin_dict[ index ] = 0
				index += resolution

			bins[ chromo ] = chr_bin_dict

	recom_file = open( recom_stats_file_name, "r" )
	skip_header = 1
	for events in recom_file: 

		if skip_header:
			skip_header = 0
			continue 

		events_list = events.strip().split(",")

		# Populate the list with the events
		event_chromosome = events_list[7]
		event_start = int( events_list[8] )
		event_length = int( events_list[9] )

		# Translates chromosomal locations into bin size
		start_bin = event_start - (event_start % resolution)
		max_bin = start_bin + (event_length - (event_length % resolution))
		number_of_bins = (max_bin - start_bin) // resolution

		prog_genotype = events_list[-1]	# Choose which of the dictionaries will be used to store the event (is the event 1, 2, or n)

		if "g" in stack_type:
			if prog_genotype == a: 
				histog_dict = dict_list[0]
			elif prog_genotype == b:  
				histog_dict = dict_list[1]
			elif prog_genotype == h: 
				histog_dict = dict_list[2]
			
			for chr_location_bin in histog_dict[event_chromosome]: 
				if (chr_location_bin >= start_bin) and (chr_location_bin < max_bin): # Add frequency of events to each chromosomal bin
					histog_dict[event_chromosome][chr_location_bin] += 1

		elif "r" in stack_type: 
			recom_bins[ event_chromosome ][ start_bin ] += 1
			recom_bins[ event_chromosome ][ max_bin ] += 1

			
	out_file = open( output_file_name, "w" )

	if "g" in stack_type: 

		for dicts in dict_list: 
			if dicts == dict_list[0]: 
				x = "Allele 1"
			elif dicts == dict_list[1]: 
				x = "Allele 2"
			elif dicts == dict_list[2]: 
				x = "Heterozygous"

			for chromosomes in dicts: 
				for bins, frequency in dicts[chromosomes].items():
					print( ",".join([x]+[str(chromosomes), str(bins), str(univ_chr_dict[chromosomes]+bins), str(frequency)]), file = out_file )

	elif "r" in stack_type: 
		
		for chromosomes in recom_bins: 
			for bins, frequency in recom_bins[chromosomes].items():
				print( ",".join([str(chromosomes), str(bins), str(univ_chr_dict[chromosomes]+bins), str(frequency)]), file = out_file )
			
	out_file.close()

def main(create_tester = 0):

	"""Main combines all the functions into one smooth workflow!"""

	arg_list = []
	for arg in sys.argv:
		arg_list.append( arg )
	
	path_and_name = arg_list[1] # Parse the path and file name

	# Split the name into two pieces: the directory path and the file name itself
	dir_path, file_name = parse_path( path_and_name )

	if (path_and_name.strip().lower() == "test"): 
		path_and_name = "test/test.csv"
		dir_path, file_name = parse_path( path_and_name )
		test_creator( path_and_name, 3, 100, 2)

	marker_cleaner( path_and_name )
	
	cleaned_file_name = dir_path + "cleaned_" + file_name
	undesirable_file_name = dir_path + "undesirables_" + file_name

	four_to_F2( cleaned_file_name )

	f2_file_name = dir_path + "f2_cleaned_" + file_name
	stats_file_name = dir_path + "statistics_" + file_name
	recombination_analysis( f2_file_name, stats_file_name )

	chr_lengths = { "Ca21chr1_C_albicans_SC5314":3190000, \
					"Ca21chr2_C_albicans_SC5314":2230000, \
					"Ca21chr3_C_albicans_SC5314":1800000, \
					"Ca21chr4_C_albicans_SC5314":1600000, \
					"Ca21chr5_C_albicans_SC5314":1190000, \
					"Ca21chr6_C_albicans_SC5314":1030000, \
					"Ca21chr7_C_albicans_SC5314":950000,  \
					"Ca21chrR_C_albicans_SC5314":2290000 }

	geno_stats_file_name = dir_path + "genos_" + file_name
	recom_Stacker( chr_lengths, stats_file_name, geno_stats_file_name, "g" ) # Stack genotypes
	os.system( "Rscript recombination_bins.R " + geno_stats_file_name + " g")

	recom_data_name = dir_path + "recom_" + file_name
	recom_Stacker( chr_lengths, stats_file_name, recom_data_name, "r" ) # Stack recombinations
	os.system( "Rscript recombination_bins.R " + recom_data_name + " r")
	

main()
