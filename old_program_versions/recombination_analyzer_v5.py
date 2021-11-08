
import statistics
import operator
import random
import numpy
import copy
import sys
import os

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

	if isinstance(csv_list, list): 
		for items in csv_list: 
			line = ",".join(items)
			print( line, file = output_file)
	elif isinstance(csv_list, dict):
		for keys, values in csv_list.items(): 
			line = ",".join([keys]+values)
			print( line, file = output_file)
	else: 
		print( "csv_printer error: This object is not a list or a dictionary" )
		quit()


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

	print_undesirables = 0
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

	if (len(chuckers)-1+len(keepers)-1) == 0:
		print( "No markers removed. Did you already prune your data? It's suspiciously clean.")
	else: 
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

	trp_f2 = transpose( f2_array )

	new_file = open( dir_path + "f2_" + file_name, "w")
	csv_printer( trp_f2, new_file )
	new_file.close()


def deploidy( f2_file_name, ploid_f2_file_name, chr_labs, blacklist_file_name = "" ): 

	""" This function removes progeny entirely and/or aneuploid chromosomes from further analysis. 
	I also don't really think that this has to be relegated to f2 data, it could definitely be put 
	anywhere in the process. """
	
	# Obtain my f2 data and the blacklist
	f2_geno_file = open( f2_file_name, "r" )
	f2_geno_array = csv_reader( f2_geno_file )
	blklist_file = open( blacklist_file_name, "r" )
	blk_list_array = csv_reader( blklist_file )
	
	# Move the data from the black list file and relabel them into a dictionary
	ploid_dict = {}
	for progeny in blk_list_array: 
		name = progeny[0].strip()
		chromosomes = progeny[1:]
		ploid_dict[ name ] = []

		# Begin adding chromosomes to the progeny
		for character in chromosomes:
			# Clean up the chromosome input from the file to make it easier to deal with, this can deffo be expanded 
			character = character.replace("\"", "")
			for long_name, abbrev_name in chr_labs.items():
				# Check to see if the label 
				if (character == abbrev_name): 
					ploid_dict[ name ].append( long_name )
					# ploid_dict[ name ].append( abbrev_name ) # Add the abbreviated name if you so choose
				elif (character == long_name): 
					ploid_dict[ name ].append( long_name )

	# In the f2 file, find indexes of markers for each chromosome (or marker numbers, though indexes will be more flexible to other pieces of data)
	chromosomes = f2_geno_array[0]
	chr_pos = f2_geno_array[1]
	mkr_number = f2_geno_array[2]
	# strains = f2_geno_array[3:] # This includes the parents... Probably shouldn't include the parents
	parent_1 = f2_geno_array[3]
	parent_2 = f2_geno_array[4]
	strains = f2_geno_array[5:] # This does not include the parents...

	chr_index_dict = {}
	for index, chr_strings in enumerate(chromosomes):
		
		# For each chromosome, get all the indexes of the markers for each chromosome
		if index == 0: 
			continue
		if chr_strings not in chr_index_dict: 
			chr_index_dict[ chr_strings ] = []
			chr_index_dict[ chr_strings ].append( index )
		else: 
			chr_index_dict[ chr_strings ].append( index )

	# Combine the dictionaries to remove the chromosomes more easily
	remove_dict = {}
	for prog_name, chr_removal_list in ploid_dict.items():
		for chromosome in chr_removal_list:
			for chr_names, indicies in chr_index_dict.items(): 
				if prog_name not in remove_dict: 
					remove_dict[ prog_name ] = chr_index_dict[ chromosome ]
				else:
					remove_dict[ prog_name ] = remove_dict[ prog_name ] + chr_index_dict[ chromosome ]

	# Construct the container for the cleaned data
	deploid_dict = {}
	deploid_dict[ chromosomes[0] ] = chromosomes[1:]
	deploid_dict[ chr_pos[0] ] = []
	for pos in chr_pos[1:]: 
		deploid_dict[ chr_pos[0] ].append( str(pos) )
	deploid_dict[ mkr_number[0] ] = []
	for num in mkr_number[1:]: 
		deploid_dict[ mkr_number[0] ].append( str(num) )
	deploid_dict[ parent_1[0] ] = parent_1[1:]
	deploid_dict[ parent_2[0] ] = parent_2[1:]

	# For each progeny, either remove it or remove the markers from the aneuploid chromosomes.
	for single_strain in strains:  
		prog_name = single_strain[0]
		
		if prog_name in ploid_dict: 
			if ploid_dict[ prog_name ] != []:
				deploid_dict[ prog_name ] = [] 
				for index, marker in enumerate(single_strain): 
					if index == 0: 
						continue
					if index in remove_dict[ prog_name ]: 
						# Add missing marker to deploid_dict[ prog_name ] if the marker is to be excluded
						deploid_dict[ prog_name ].append( "-" )
					else: 
						deploid_dict[ prog_name ].append( str(marker) )
			else: 
				continue # Remove the progeny. Continue is not necessary, but if anybody wants to add something new, this is the place. 
		else: 
			deploid_dict[ prog_name ] = single_strain[1:]

	ploid_f2_file = open( ploid_f2_file_name, "w" )
	csv_printer( deploid_dict, ploid_f2_file )
	ploid_f2_file.close()

	number_o_prog = len( deploid_dict )-5 

	return deploid_dict, number_o_prog


def marker_tally_ho( strain_marker_dictionary, stack_bar_file_name, a = "1", b = "2", h = "n", missing = "-", include_missing = 0 ):

	""" Tallys the markers for each strain and outputs a file for creating a stacked bar chart. """ 

	stack_bar_file = open( stack_bar_file_name, "w" )

	if include_missing: 
		header = "Progeny/Strain, allele, count, a_prop, b_prop, h_prop, m_prop"
		print( header, file = stack_bar_file )

	else: 
		header = "Progeny/Strain, allele, count, a_prop, b_prop, h_prop"
		print( header, file = stack_bar_file )

	strain_tally_dict = {}
	for strain_names, marker_list in strain_marker_dictionary.items(): 

		if  (strain_names == "Chromosome") or \
			(strain_names == "Chr_Position") or \
			(strain_names == "# Catalog ID"):
			continue
		strain_tally_dict[ strain_names ] = { a:0, b:0, h:0, missing:0 } 
		strain_tals = strain_tally_dict[ strain_names ]

		for markers in marker_list: 

			if markers == a: 
				strain_tals[ a ] += 1
			elif markers == b:
				strain_tals[ b ] += 1
			elif markers == h:
				strain_tals[ h ] += 1
			elif markers == missing:
				strain_tals[ missing ] += 1
			else: 
				print( "Unrecognized marker in ", strain_names )
				print( "Quitting" )
				quit()

		if include_missing:
			total = strain_tals[a] + strain_tals[b] + strain_tals[h] + strain_tals[missing]
			
			a_proportion = round( strain_tals[a]/total, 2 )
			b_proportion = round( strain_tals[b]/total, 2 )
			h_proportion = round( strain_tals[h]/total, 2 )
			m_proportion = round( strain_tals[missing]/total, 2 )

			for alleles, count in strain_tals.items(): 
				line = ",".join([ str(strain_names), str(alleles), str(count), str(a_proportion), str(b_proportion), str(h_proportion), str(m_propotion) ])
				print( line, file = stack_bar_file )

		else: 
			total = strain_tals[a] + strain_tals[b] + strain_tals[h]

			a_proportion = round( strain_tals[a]/total, 2 )
			b_proportion = round( strain_tals[b]/total, 2 )
			h_proportion = round( strain_tals[h]/total, 2 )

			for alleles, count in strain_tals.items(): 
				if alleles == missing: 
					continue
				line = ",".join([ str(strain_names), str(alleles), str(count), str(a_proportion), str(b_proportion), str(h_proportion) ])
				print( line, file = stack_bar_file )

	stack_bar_file.close()


def recombination_analysis( f2_file_name, output_file_name, centromeres, chr_lengths, a = "1", b = "2", h = "n", centromere_rules = 1):
	
	"""This function takes the cleaned f2 data and searches for recombination events, and stores them. """

	f2_file = open( f2_file_name, "r" )
	output_file = open( output_file_name, "w" )

	f2_array = csv_reader( f2_file )

	# Printing the stats file
	header_str = "Sample Name, Recombined Chromosome, Recombination Start Site (bp), Recombination Event Length (bp)"
	print( header_str , file = output_file )

	chromosomes = f2_array[0]
	chr_pos     = f2_array[1]
	mrkr_cat_id = f2_array[2]
	parent1     = f2_array[3]
	parent2     = f2_array[4]
	children    = f2_array[5:]
	
	progs_dict = {}
	recombs = {}
	left_rite_indies = {}
	recom_count = 0

	# Split the fucking chromosomes at the centromeres... 
	for chro, centros in centromeres.items(): 
		left_arm_key = chro + "_left"
		rite_arm_key = chro + "_rite"
		left_rite_indies[ left_arm_key ] = []
		left_rite_indies[ rite_arm_key ] = []
		centromere_left_bound = centros[0]
		centromere_rite_bound = centros[1]
		center = (centromere_left_bound + centromere_rite_bound)/2
		for i, chr_strings in enumerate(chromosomes):
			if chro == chr_strings: 
				if int(chr_pos[i]) <= center:
					left_rite_indies[ left_arm_key ].insert( 0,i ) # These will be listed to count recombinations away from the centromere
				if int(chr_pos[i]) >= center:
					left_rite_indies[ rite_arm_key ].append( i ) # In order they appear in the 4 way file

	for progeny in children: 
		info_sub_dict = {}
		for chr_arm, index_list in left_rite_indies.items(): 
			info_sub_dict[chr_arm] = { "marker": [], "chr_pos": [] }
			for index in index_list: 
				info_sub_dict[chr_arm]["marker" ].append( progeny[index] )
				info_sub_dict[chr_arm]["chr_pos"].append( int(chr_pos[index]) )
		progs_dict[progeny[0]] = info_sub_dict


	# Begin assessing recombination:
	for prog_name, chr_arm_dict in progs_dict.items(): 

		recombs[ prog_name ] = { "marker_tally":{}, "recoms": [] }
		prog_data = recombs[ prog_name ]

		for chr_arm, mrkr_and_pos in chr_arm_dict.items(): 
			chr_name = chr_arm[:-5]
			marker_list = mrkr_and_pos[ "marker" ]
			posit_list = mrkr_and_pos[ "chr_pos" ]

			stor_idx = "empty"
			first_recom = 1
			for index in range(0,len(marker_list)):

				marker = marker_list[ index ]

				prog_data[ "marker_tally" ][ ",".join([chr_name, str(posit_list[index])])] = marker
				if progeny[i] == "-":
					continue

				if stor_idx == "empty": # First time around
					stor_idx = index  

				if first_recom == 1 and (marker == marker_list[stor_idx]):
					stor_idx = index # Keep moving the stored index until the first time something is different. 
				
				elif first_recom == 1 and (marker != marker_list[stor_idx]): 
					first_recom = 0
					recom_count += 1
					# evaluate size and print to a new file
					event_size = abs(int(posit_list[stor_idx]) - posit_list[index])
					event = ",".join( [chr_name, str(posit_list[stor_idx]), str(event_size)] )
					print( prog_name + "," + event, file = output_file)
					prog_data[ "recoms" ].append( event )				
					stor_idx = index 

				elif first_recom == 0 and (marker != marker_list[stor_idx]):

					recom_count += 1
					event_size = abs(int(posit_list[stor_idx]) - posit_list[index])
					event = ",".join( [chr_name, str(posit_list[stor_idx]), str(event_size)] )
					print( prog_name + "," + event, file = output_file)
					prog_data[ "recoms" ].append( event )
					stor_idx = index

				elif first_recom == 1 and ((index + 1) == len(marker_list)):
					pass # No recombination

				elif first_recom == 0 and ((index + 1) == len(marker_list)):

					recom_count += 1
					if chr_arm[-4:] == "left": 
						event_size = abs(int(posit_list[stor_idx]))
					elif chr_arm[-4:] == "rite": 
						event_size = abs(int(posit_list[stor_idx]) - chr_lengths[chr_name])
					event = ",".join( [chr_name, str(posit_list[stor_idx]), str(event_size)] )
					print( prog_name + "," + event, file = output_file)
					prog_data[ "recoms" ].append( event )
					stor_idx = index
	
	print( recombs )
	print( recom_count )
	return recombs # Recombs is just a dictionary with progeny names as keys and then a list of strings with recombination "events" 
	# Sample event: event = ",".join( [chromosomes[i], str( left_side_centro ), str(event_size)] )

	f2_file.close()
	output_file.close()



def recom_Stacker( chr_lengths, recom_dict, output_file_name, resolution = 100000, a = "1", b = "2", h = "n", missing = "-" ): 
	
	"""Creates a dictionary of chromosomes and bin sizes that wll be used for histograms to depict recombination hotspot for each allele type."""   

	# Compute the universal position for each chromosome into the total variable. 
	# This is to translate marker positions into universal positions for graphing 
	total = 0
	univ_chr_dict = {}
	for chromosomes, lengths in chr_lengths.items(): 
		univ_chr_dict[ chromosomes ] = total
		total += int(lengths)

	# Build a dictionary to store the recombination events and for tallying the markers
	recom_bins = {} 
	tally_dict = {} 
	for chromo, length in chr_lengths.items(): 
		chr_bin_dict = {}
		tal_bin_dict = {}
		for i in range(0, length-(length%resolution)+resolution, resolution):
			chr_bin_dict[ i ] = 0
			tal_bin_dict[ i ] = { a: 0, b: 0, h: 0, missing: 0 }
		recom_bins[ chromo ] = chr_bin_dict
		tally_dict[ chromo ] = tal_bin_dict

	# I made two loops instead of one because handling the data from both simultaneously would be extraordinarily difficult.  
	for progeny_name, marker_and_recom_dict in recom_dict.items(): 

		# First of two loops; this one for recombination events
		for chr_pos_l1 in marker_and_recom_dict[ "recoms" ]:

			# "recoms" has chromosome,start,end list format for each individual progeny
			cpl1 = chr_pos_l1.split(",")
			chromosome, position, event_size = cpl1[0], int(cpl1[1]), int(cpl1[2]),
			
			# This line will drop the position into a given bin by removing the remainder of resolution
			# start_bin = position - (position % resolution)  
			end_bin = position + event_size - ( (position + event_size) % resolution )
			
			# Put that into the appropriate bin of tally_dict:
			# recom_bins[ chromosome ][ start_bin ] += 1
			recom_bins[ chromosome ][ end_bin ] += 1

		# Second loop for the marker tallying
		for chr_pos_l2, marker in marker_and_recom_dict[ "marker_tally" ].items(): 
			
			# Marker_tally has chromosome,xxxxx key format for each individual progeny
			chromosome, position = chr_pos_l2.split(",")[0], int(chr_pos_l2.split(",")[1]) 
			
			# With the position data, compute which bin the marker should go into on the regular gene coordinates
			t_bin = position - (position % resolution)
			
			# Put that into the appropriate bin of tally_dict: 
			tally_dict[ chromosome ][ t_bin ][ marker ] += 1

	# Print the information
	recom_out_file = open( output_file_name, "w" )
	header = "chr,chr.pos,univ.pos,a,b,h,missing,recoms"  
	print( header, file = recom_out_file )
	for chromosome, chr_bins in recom_bins.items(): 
		for bins, recom_freq in chr_bins.items(): 
			# I'm gonna call the tallies using the recombination data because the positions are identical and I want them in the same file. 
			prt_pos = str(bins)
			prt_univ_pos = str(bins + univ_chr_dict[ chromosome ])
			prt_a = str(tally_dict[ chromosome ][ bins ][ a ])
			prt_b = str(tally_dict[ chromosome ][ bins ][ b ])
			prt_h = str(tally_dict[ chromosome ][ bins ][ h ])
			prt_m = str(tally_dict[ chromosome ][ bins ][ missing ])
			prt_freq = str(recom_freq)

			print( ",".join( [ chromosome, prt_pos, prt_univ_pos, prt_a, prt_b, prt_h, prt_m, prt_freq ] ), file = recom_out_file )

	recom_out_file.close()

	return recom_bins, resolution 


def centromerely_math( chr_lengths, centromere_locations, recom_dict, output_file_name, \
						resolution = 100000, count_incentro_bool = 0, left_arm = "left", right_arm = "right", centromere = "centromere" ): 

	""" Calculate cumulative recombination as a function of distance from centromeres """

	# Compute the universal position for each chromosome into the total variable. 
	# This is to translate marker positions into universal positions for graphing 
	total = 0
	univ_chr_dict = {}
	for chromosomes, lengths in chr_lengths.items(): 
		univ_chr_dict[ chromosomes ] = total
		total += int(lengths)

	# Make a dictionary of positions with recombination events. 
	centro_recoms = {} 
	sorted_recoms = {}
	cyuum_recoms = {}
	for progeny, mkr_data in recom_dict.items(): 
		for recoms in mkr_data["recoms"]: 
			recom_list = recoms.split(",")
			chro, recom_start, recom_size = recom_list[0], int(recom_list[1]), int(recom_list[2])
			recom_end = recom_start + recom_size
			
			if chro not in centro_recoms:  # Using recom_end because that's what I use to count recombination for each position
				centro_recoms[ chro ] = { left_arm:{}, right_arm:{}, centromere:{} }
				sorted_recoms[ chro ] = { left_arm:[], right_arm:[], centromere:[] } # This dictionary and the next are for later
				cyuum_recoms[ chro ] = { left_arm:{}, right_arm:{}, centromere:0 }


			# Make sure the recombination isn't anywhere else
			put_er_in_bool = (recom_end not in centro_recoms[ chro ][ left_arm   ]) and \
			                 (recom_end not in centro_recoms[ chro ][ right_arm ]) and \
			                 (recom_end not in centro_recoms[ chro ][ centromere  ])

			# If the recombination can be put into the dictionary, then put it in there. 
			if put_er_in_bool: 
				if recom_end <= centromere_locations[ chro ][0]: 
					centro_recoms[ chro ][right_arm][ recom_end ] = 1
				elif recom_end >= centromere_locations[ chro ][1]: 
					centro_recoms[ chro ][left_arm][ recom_end ] = 1
				else: 
					centro_recoms[ chro ][ centromere ][ recom_end ] = 1

			# If it's already in there, then add another event to it
			else:
				if recom_end <= centromere_locations[ chro ][0]: 
					centro_recoms[ chro ][right_arm][ recom_end ] += 1
				elif recom_end >= centromere_locations[ chro ][1]: 
					centro_recoms[ chro ][left_arm][ recom_end ] += 1
				else: 
					centro_recoms[ chro ][ centromere ][ recom_end ] += 1

	# All this just to sort the chromosomal positions... I should add a centromere point as zero....
	for chromosomes, centro_rel_pos in centro_recoms.items():
		for chrom_arm, positions in centro_rel_pos.items():
			empty_list = []
			centro_loc = round(statistics.mean(centromere_locations[ chromosomes ]))
			empty_list.append( [ centro_loc, 0 ])
			
			for position, recom_count in positions.items(): 
				empty_list.append( [ position, recom_count ] )

			# No need to sort the centromere stuff, if anything is inside it because those are treated differently
			if chrom_arm == left_arm:
				sorted_recoms[ chromosomes ][ chrom_arm ] = sorted(empty_list, key = operator.itemgetter(0))
			elif chrom_arm == right_arm:
				sorted_recoms[ chromosomes ][ chrom_arm ] = sorted(empty_list, key = operator.itemgetter(0), reverse = True)
	
	# With the sorted positions, I can now quickly go through the data and find the cumulative values as distance from the chromosome increases
	for chromosomes, centro_rel_pos in sorted_recoms.items():
		
		# Establish the number of recombination events in the centromere for reference at all other locations
		centro_total = 0
		if count_incentro_bool:
		# in-centromere recombination events will be counted (not common, and/or impossible), so kept out by default
			for recoms in centro_rel_pos[ centromere ]:
				centro_total += recoms[1]
				cyuum_recoms[ chromosomes ][ centromere ] = centro_total
		else: 
			cyuum_recoms[ chromosomes ][ centromere ] = centro_total  

		# Accumulate recombination events from the lowest chromosomal position to the highest
		upstream_total = centro_total 
		for recoms in centro_rel_pos[ left_arm ]:
			cyuum_recoms[ chromosomes ][ left_arm ][ recoms[0] ] = upstream_total + recoms[1]
			upstream_total += recoms[1]
			 
		# Accumulate recombination events from the highest to the lowest
		downstream_total = centro_total 	
		for recoms in centro_rel_pos[ right_arm ]:
			cyuum_recoms[ chromosomes ][ right_arm ][ recoms[0] ] = downstream_total + recoms[1]
			downstream_total += recoms[1]
	
	# Put all the events together to be written into recombination_bins.R
	final_dict = {}
	for chromosomes, centro_rel_pos in cyuum_recoms.items():
		final_dict[ chromosomes ] = []
		# Centromere counts are factored into the upstream and downstream counts
		for recom_pos, recom_count in centro_rel_pos[ right_arm ].items():
			final_dict[ chromosomes ].append( [ recom_pos, recom_count, left_arm ] )
		for recom_pos, recom_count in centro_rel_pos[ left_arm ].items():
			final_dict[ chromosomes ].append( [ recom_pos, recom_count, right_arm ] ) 

	centro_file = open( output_file_name, "w" )
	header = "chromosome,position,absolute_pos,cumulative_recom_cnt,chromosome_arm "
	print( header, file = centro_file )
	for chromosomes in final_dict: 
		final_dict[ chromosomes ] = sorted( final_dict[ chromosomes ], key = operator.itemgetter(0) )

		for events in final_dict[ chromosomes ]: 
			absolute_pos = events[ 0 ] + univ_chr_dict[ chromosomes ]
			accumulated_recom_count = events[1]
			arm = events[2]
			print( ",".join([ chromosomes, str(events[0]), str(absolute_pos), str(accumulated_recom_count), events[2] ]), file = centro_file )


def main(create_tester = 0):

	"""Main combines all the functions into one "smooth" workflow!"""

	arg_list = []
	for arg in sys.argv:
		arg_list.append( arg )
	
	path_and_name = arg_list[1] # Parse the path and file name

	# Split the name into two pieces: the directory path and the file name itself
	dir_path, file_name = parse_path( path_and_name )

	# Taken from the genes.gff from A21 
	centromere_locations = { "Ca21chr1_C_albicans_SC5314": [1563038,1565967], \
							 "Ca21chr2_C_albicans_SC5314": [1927255,1930214], \
							 "Ca21chr3_C_albicans_SC5314": [823333,826481], \
							 "Ca21chr4_C_albicans_SC5314": [992579,996216], \
							 "Ca21chr5_C_albicans_SC5314": [468716,471745], \
							 "Ca21chr6_C_albicans_SC5314": [980040,983792], \
							 "Ca21chr7_C_albicans_SC5314": [425812,428712], \
							 "Ca21chrR_C_albicans_SC5314": [1743190,1747664] }
	
	# Size of each chromosome in basepairs. Exact is best, but these can be estimated. 
	chr_lengths = { "Ca21chr1_C_albicans_SC5314":3190000, \
					"Ca21chr2_C_albicans_SC5314":2230000, \
					"Ca21chr3_C_albicans_SC5314":1800000, \
					"Ca21chr4_C_albicans_SC5314":1600000, \
					"Ca21chr5_C_albicans_SC5314":1190000, \
					"Ca21chr6_C_albicans_SC5314":1030000, \
					"Ca21chr7_C_albicans_SC5314":950000,  \
					"Ca21chrR_C_albicans_SC5314":2290000 }

	chr_labs =  {	"Ca21chr1_C_albicans_SC5314": "1", \
					"Ca21chr2_C_albicans_SC5314": "2", \
					"Ca21chr3_C_albicans_SC5314": "3", \
					"Ca21chr4_C_albicans_SC5314": "4", \
					"Ca21chr5_C_albicans_SC5314": "5", \
					"Ca21chr6_C_albicans_SC5314": "6", \
					"Ca21chr7_C_albicans_SC5314": "7", \
					"Ca21chrR_C_albicans_SC5314": "R" }

	# If performing a test run, this will be the section that is activated
	test_bool = path_and_name.strip().lower() == "test"
	if test_bool: 
		path_and_name = "test/test.csv"
		dir_path, file_name = parse_path( path_and_name )
		test_creator( path_and_name, 1, 100, 500)

	marker_cleaner( path_and_name )
	
	cleaned_file_name = dir_path + "cleaned_" + file_name
	undesirable_file_name = dir_path + "undesirables_" + file_name

	four_to_F2( cleaned_file_name )

	f2_file_name = dir_path + "f2_cleaned_" + file_name
	ploid_f2_file_name = dir_path + "f2_recleaned_" + file_name
	if test_bool:
		blacklist_file_name = "test/blk_list_test.csv"
	else:
		blacklist_file_name = "progeny_ploidy_blacklist.csv"
		# blacklist_file_name = "blank_blacklist.csv" # Use this line to compare "before and after" blacklisting
	
	deploid_dict, number_o_progeny = deploidy( f2_file_name, ploid_f2_file_name, chr_labs, blacklist_file_name )
	stack_bar_file_name = dir_path + "stackbar_" + file_name 
	marker_tally_ho( deploid_dict, stack_bar_file_name )
	os.system( "Rscript recombination_bins.R " + stack_bar_file_name + " s" )

	stats_file_name = dir_path + "statistics_" + file_name
	recom_dict = recombination_analysis( ploid_f2_file_name, stats_file_name, centromere_locations, chr_lengths )
	os.system( "Rscript recombination_bins.R " + stats_file_name + " h" )

	resolve_this = 100000
	recom_data_name = dir_path + "recom_" + file_name
	recom_bins, resolution = recom_Stacker( chr_lengths, recom_dict, recom_data_name, resolution = resolve_this ) # Stack recombinations

	os.system( "Rscript recombination_bins.R " + recom_data_name + " g")

	centro_data_name = dir_path + "centro_" + file_name
	centromerely_math( chr_lengths, centromere_locations, recom_dict, centro_data_name, resolution = resolve_this ) # takes The dictionary from recombination_analysis(), not recom_Stacker()

	os.system( "Rscript recombination_bins.R " + centro_data_name + " c" )

main()
