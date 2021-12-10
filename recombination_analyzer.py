
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


def marker_cleaner(dir_path, file_name, blacklist_file_name = "remove_markers.csv", \
					remove_strangers = 1, blacklisted = 1, print_undesirables = 1): 
	"""This function removes markers that are missing from eiter parent and markers that are
	ambiguous between parents (the parents share an allele). """

	four_way_file_name = dir_path + file_name

	four_way_file = open( four_way_file_name, "r" )
	untrp_leels = csv_reader( four_way_file ) # Load the file, but it needs to be transposed
	leels = transpose( untrp_leels )

	try: 
		if blacklisted:
			blacklist_file = open( dir_path + blacklist_file_name, "r" )
			blk_mkr_list = []
			for marker_cat_ID in blacklist_file: 
				if "#" in marker_cat_ID: 
					continue
				blk_mkr_list.append( marker_cat_ID.strip() )
			blacklist_file.close()

	except IOError: 
		blk_mkr_list = []
		print( "No black-listed marker file; making one: ", dir_path + blacklist_file_name )
		blacklist_file = open( dir_path + blacklist_file_name, "w" )
		print( "# This file is for removing markers that STACKS miscalled; use IGV to manually curate \"suspicious\"",\
		" (highly recombinant) markers from f2_recleaned_" + file_name + ".", file = blacklist_file)

	if blacklisted: 
		before = len(leels) 
		for markers in leels: 
			marker_num = markers[2].strip()
			if marker_num in blk_mkr_list: 
				leels.remove(markers)
		after =  len(leels) 
		blk_lst_rmv = before - after

	keepers = []
	chuckers = []
	strangers = []

	stranger_count = 0

	header = 1
	for lines in leels: 
		if header: 
			keepers.append( lines )
			chuckers.append( lines )
			strangers.append( lines )
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

		unique_alleles = set()

		unique_alleles.add(parent1[0])
		unique_alleles.add(parent1[1])
		unique_alleles.add(parent2[0])
		unique_alleles.add(parent2[1])

		parental_leel_count = len( unique_alleles )

		if (parental_leel_count > 2) and ( remove_strangers ): # Remove tandem repeat recombination markers
			strangers.append( lines )
			stranger_count += 1 
			continue

		keepers.append( lines )

	if (len(chuckers)-1+len(keepers)-1) == 0:
		print( "No markers removed. Did you already prune your data? It's suspiciously clean.")
	else: 
		if blacklisted:
			print( "marker_cleaner removed ", len(chuckers) - 1 + stranger_count , \
			# " markers (", 100*(len(chuckers)-1+ stranger_count)/(len(chuckers)-1+len(keepers)-1+stranger_count), "%","removed )\n", \
			# "in addition to", blk_lst_rmv, "black-listed markers." )
			" markers (", str(100*(len(chuckers)-1+ stranger_count)/(len(chuckers)-1+len(keepers)-1+stranger_count)) + "%" + " removed )" )
		else:
			print( "marker_cleaner removed ", len(chuckers) - 1 + stranger_count, \
				" markers (", str(100*(len(chuckers)-1+ stranger_count)/(len(chuckers)-1+len(keepers)-1+stranger_count)) + "%" + " removed )")

	srt_keepers = sorted( keepers, key=operator.itemgetter(0,1) )
	for lines in srt_keepers: 
		lines[1] = str(lines[1]) # Convert back to strings so they can be output. 
		lines[2] = str(lines[2])
	
	header = srt_keepers.pop(-1) # Move the header string back to the front
	srt_keepers.insert(0, header)

	# Remove black-listed markers here.
	trp_keepers = transpose( srt_keepers ) # File needs to be tranposed again. 
	new_file = open( dir_path + "cleaned_" + file_name, "w" )
	for lines in trp_keepers: 
		print( ",".join(lines), file = new_file )
	new_file.close()

	if print_undesirables:
		srt_strangers = sorted( strangers, key = operator.itemgetter(0,1) )
		srt_chuckers = sorted( chuckers, key=operator.itemgetter(0,1) )
		for lines in srt_chuckers:
			lines[1] = str(lines[1]) # Convert back to strings so they can be output. 
			lines[2] = str(lines[2]) 

		for lines in srt_strangers:
			lines[1] = str(lines[1])
			lines[2] = str(lines[2])
		
		header = srt_chuckers.pop(-1) # Move the header string back to the front
		stranger_header = srt_strangers.pop(-1)
		srt_chuckers.insert(0, header)
		srt_strangers.insert(0, stranger_header)

		trp_chuckers = transpose( srt_chuckers )
		trp_strangers = transpose( srt_strangers )
		new_str_file = open( dir_path + "strangers_" + file_name, "w" ) 
		new_un_file = open( dir_path + "undesirables_" + file_name, "w" )

		for lines in trp_chuckers: 
			print( ",".join(lines), file = new_un_file )
		new_un_file.close()
		for lines in trp_strangers: 
			print( ",".join(lines), file = new_str_file )
		new_str_file.close()


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


def deploidy( dir_path, f2_file_name, ploid_f2_file_name, chr_labs, blacklist_file_name = "", \
				mark_strangers = 1, replace_strange = "-", quit_on_no_strange = 0 ): 

	""" This function removes entire progeny, aneuploid chromosomes, or individual markers (after manually checking them) from further analysis. 
	I also don't really think that this has to be relegated to f2 data, it could definitely be put 
	anywhere in the process. """
	
	# Obtain my f2 data and the blacklist
	f2_geno_file = open( dir_path + f2_file_name, "r" )
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
				for index, marker in enumerate(
					single_strain): 
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

	try:
		strange_marker_file = open( dir_path + "strange_markers.csv", "r" )
		strange_mkr_list = csv_reader( strange_marker_file )
		strange_cat_nums = strange_mkr_list[2]
	
	except IOError: 
		if quit_on_no_strange: 
			quit()
		else: 
			print("Could not find " + dir_path + "strange_markers.csv. We wrote an empty file to fill." )
			# strange_marker_file = open( dir_path + "strange_markers.csv", "w" )
		mark_strangers = 0 

	except IndexError: 
		print( "No manually curated markers in " + dir_path + "strange_markers.csv. Default behavior is to continue without them.")
		mark_strangers = 0 
	
	if mark_strangers:	

		strange_prog_dict = {}
		for csv_lines in strange_mkr_list: 
			prog_name = csv_lines[0]
			strange_prog_dict[ prog_name ] = csv_lines[1:]
		# Find indexes in cleaned_markers where strange markers are
		strange_idx_list = []
		for index in range(1,len(mkr_number)): 
			if int(mkr_number[ index ]) in strange_cat_nums:
				strange_idx_list.append(index-1)

		for idx in strange_idx_list: 
			for prog_name, mkr_list in deploid_dict.items():

				# Parental names shouldn't be too much of a problem because they shouldn't contain any x's
				if  ("529L_" in prog_name)        or \
					("SC5314" in prog_name)       or \
					("Chromosome" in prog_name)   or \
					("Chr_Position" in prog_name) or \
					("# Catalog ID" in prog_name): 
					continue

				new_marker_index = strange_idx_list.index(idx)
				new_marker = strange_prog_dict[prog_name][new_marker_index].strip()

				if mkr_list[ idx ] == "-": 
					continue
				else: 
					if new_marker == "x": 
						mkr_list[ idx ] = replace_strange 
					else: 
						mkr_list[ idx ] = new_marker 

	ploid_f2_file = open( dir_path + ploid_f2_file_name, "w" )
	csv_printer( deploid_dict, ploid_f2_file )
	ploid_f2_file.close()

	number_o_prog = len( deploid_dict )-5 

	return deploid_dict, number_o_prog


def marker_tally_ho( strain_marker_dictionary, stack_bar_file_name, \
						a = "1", b = "2", h = "n", missing = "-", strange = "x", include_missing = 0 ):

	""" Tallys the markers for each strain and outputs a file for creating a stacked bar chart. """ 

	stack_bar_file = open( stack_bar_file_name, "w" )
 
	header = "Progeny/Strain, allele, count, prop, a_prop, b_prop, h_prop"
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
			markers = markers.strip()
			if markers == a: 
				strain_tals[ a ] += 1
			elif markers == b:
				strain_tals[ b ] += 1
			elif markers == h:
				strain_tals[ h ] += 1
			elif markers == strange: 
				strain_tals[ h ] += 1 # Technically, strange markers are heterozygous markers
			elif markers == missing:
				strain_tals[ missing ] += 1
			else: 
				print( "Unrecognized marker: ", markers, " in ", strain_names )
				print( "Quitting" )
				quit()

		total = strain_tals[a] + strain_tals[b] + strain_tals[h] 
		# total = strain_tals[a] + strain_tals[b] + strain_tals[h] + strain_tals[missing]

		a_proportion = round( strain_tals[a]/total, 5 )
		b_proportion = round( strain_tals[b]/total, 5 )
		h_proportion = round( strain_tals[h]/total, 5 )
		# m_proportion = round( strain_tals[missing]/total, 5 )

		for alleles, count in strain_tals.items(): 
			if (include_missing == 0) and (missing in alleles):
				continue 
			leel_prop = round( strain_tals[alleles]/total, 5 )
			line = ",".join([ str(strain_names), str(alleles), str(count), str(leel_prop), str(a_proportion), str(b_proportion), str(h_proportion) ])
			# line = ",".join([ str(strain_names), str(alleles), str(count), str(leel_prop), str(a_proportion), str(b_proportion), str(h_proportion), str(m_proportion) ])
			print( line, file = stack_bar_file )

	stack_bar_file.close()


def recombination_analysis( dir_path, f2_file_name, output_file_name, count_filename, centromeres, chr_lengths, \
							a = "1", b = "2", h = "n", missing = "-"):
	
	"""This function takes the cleaned f2 data (data with only markers with and searches for recombination events, and stores them. """

	f2_file = open( dir_path + f2_file_name, "r" )
	f2_array = csv_reader( f2_file )

	count_file = open( count_filename, "w" )
	header = "prog.name, recombinations"
	print( header, file = count_file )

	# Printing the stats file
	output_file = open( output_file_name, "w" )
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

	# Split the chromosomes at the centromeres into left and right (rite) arms
	for chro, centros in centromeres.items(): 
		left_arm_key = chro + "_left"
		rite_arm_key = chro + "_rite"
		left_rite_indies[ left_arm_key ] = []
		left_rite_indies[ rite_arm_key ] = []
		centromere_left_bound = centros[0]
		centromere_rite_bound = centros[1]
		center = (centromere_left_bound + centromere_rite_bound)/2
		centros.append( center )

		for i, chr_strings in enumerate(chromosomes):
			if chro == chr_strings: 
				if int(chr_pos[i]) <= center:
					left_rite_indies[ left_arm_key ].insert( 0,i ) # These will be listed to count recombinations away from the centromere
				if int(chr_pos[i]) >= center:
					left_rite_indies[ rite_arm_key ].append( i ) # In order they appear in the 4 way file

	# Make a data structure to house recombination information each progeny will have chromosome arms 
	# that contain marker genotype ("marker") and position ("chr_pos") information from the parents to be used for recombination analysis. 
	for progeny in children: 
		info_sub_dict = {}
		for chr_arm, index_list in left_rite_indies.items(): 

			info_sub_dict[chr_arm] = { "marker": [], "chr_pos": [] }
			for index in index_list: 

				try:
					info_sub_dict[chr_arm]["marker" ].append( progeny[index] ) # Marker genotype identity
					info_sub_dict[chr_arm]["chr_pos"].append( int(chr_pos[index]) ) # Turning them into integers allows booleans for control flow
				except IndexError: 
					print( "Index: ", index )
					print( "Progeny length: ", len(progeny))
					print( "Progeny:", progeny )

		progs_dict[progeny[0]] = info_sub_dict

	# Begin assessing recombination:
	bir = 0 		  # bir = break-induced replication: recombinations that go to the end of the chromosomes
	recom_count = 0   # The total number of recombinations across all progeny
	total_progeny = 0 # The total number of progeny

	for prog_name, chr_arm_dict in progs_dict.items(): 

		# Remove individual missing markers (and respecitive locations) from progeny. Without this, the algorithm calls false recombinations. 
		for chr_arms, mkr_n_pos_dict in progs_dict[prog_name].items(): 
			new_mkr_list = []
			new_pos_list = []
			for index in range(0,len(mkr_n_pos_dict["marker"])): 
				if missing not in mkr_n_pos_dict["marker"][index]: 
					new_mkr_list.append(mkr_n_pos_dict["marker"][index])
					new_pos_list.append(mkr_n_pos_dict["chr_pos"][index])
		mkr_n_pos_dict["marker"] = new_mkr_list
		mkr_n_pos_dict["chr_pos"] = new_pos_list

		prog_recom_count = 0 # Total recombinations in each progeny

		total_progeny += 1 # Total number of progeny 

		# Initiate the nested dictionary for tallying markers and recording recombination events
		recombs[ prog_name ] = { "marker_tally":{}, "recoms": [] } 

		# Setting recombs[ prog_name ] as a variable for easier reading and retrieval
		prog_data = recombs[ prog_name ]  

		for chr_arm, mrkr_and_pos in chr_arm_dict.items(): 

			# Making variables for chromosome arms for variability
			chr_name = chr_arm[:-5]
			arm_side = chr_arm[-4:]

			# Making variables for markers and chromosome position for readability 
			marker_list = mrkr_and_pos[ "marker" ]
			posit_list = mrkr_and_pos[ "chr_pos" ]

			# Initiate sub-dictionary for recombination events
			recom_frag_dict = {"positions": [], "markers": []}

			# Initiate reference index variable for comparison for capturing recombination events
			stor_idx = "empty"

			# Triggers a different response after the first recombination event is found
			first_recom = 1
			
			# Begin looking for recombination events 
			for index in range(0,len(marker_list)):

				# Boolean checking for the last marker on the chromosome arm; used for special cases
				last_one = (index + 1)  == len(marker_list)

				# Make a variable for the marker; if this changes from stor_mkr, a recombination has occurred
				marker = marker_list[ index ]

				# prog_data[ "marker_tally" ] is used to make lines for an output file that will be used for allele contribution of progeny	
				prog_data[ "marker_tally" ][ ",".join([chr_name, str(posit_list[index])])] = marker

				if stor_idx == "empty": # First time around
					stor_idx = index  

				try: 
					recom_bool = (marker != marker_list[stor_idx]) and (marker != missing) # (marker != missing) is unnecessary, but this line was written before I started just removing missing markers and I don't want to break my code. It's working fine with it in here.
				except IndexError: # stor_idx is too big for marker_list 
					print( "Disjunction between variable \"stor_idx\" and the index on a given chromosome. This code can't handle your input." ) 
					recom_bool = 0
				except TypeError: # No recombinations yet, so stor_idx hasn't been changed yet
					recom_bool = 0

				# recom_frag_dict records positions and markers of a recombination fragment
				recom_frag_dict[ "positions" ].append( posit_list[index] )
				recom_frag_dict[ "markers" ].append( marker )

				if last_one: # If it's the last item in the list; this is a special case. Regular recombiantions are below
					
					if first_recom: # The first change in the marker happens at the last possible position; chromosome only changed once at the very end
						if recom_bool:
							# Tally the recombination events
							recom_count += 1
							prog_recom_count += 1

							# The main difference between left and right arms is that the chromosome coordinates run backwards in the left arm
							# Thus, the different arms require different ways to measure the recombiantion tract for the last recombinations
							if arm_side == "left":
								# Left arm starts at zero; the last marker's position is also its distance from the end of the chromosome
								event_size = abs( recom_frag_dict["positions"][-1] )
							else: 
								# Right arm ends at a non-zero length; distance to end needs to be subtracted from the total size. 
								event_size = abs( recom_frag_dict["positions"][-1] - chr_lengths[ chr_name ])
							
							# Store the event information for use of analysis: 
								# chr_name:                              name of the chromosome
								# str(recom_frag_dict["positions"][-1]): last recorded position that triggered recom_bool
								# event_size:                            distance from the change to the (for this case) end of the chromosome
							event = ",".join( [chr_name, str(recom_frag_dict["positions"][-1]), str(event_size)] )
							print( prog_name + "," + event, file = output_file)
							prog_data[ "recoms" ].append( event )

							# Increment break-induced recombination events
							bir += 1
					
					else: # There have already been recombinations along the chromosome; special case for end-of-chromosome recombinations
						if recom_bool:
							# There are two diferent recombinations happening; one before the last marker and one AT the last marker
							recom_count += 2
							prog_recom_count += 2

							# First recombination is measured from the second-to-last marker to the last marker
							event_size = abs(recom_frag_dict["positions"][0] - recom_frag_dict["positions"][-1])
							event = ",".join( [chr_name, str(recom_frag_dict["positions"][0]), str(event_size)] )
							print( prog_name + "," + event, file = output_file)
							prog_data[ "recoms" ].append( event )

							# Second recombiantion is measured from the last marker to the end of the chromosome
							if arm_side == "left":
								event_size = abs( recom_frag_dict["positions"][-1] )
							else: 
								event_size = abs( recom_frag_dict["positions"][-1] - chr_lengths[ chr_name ])
							event = ",".join( [chr_name, str(recom_frag_dict["positions"][-1]), str(event_size)] )
							print( prog_name + "," + event, file = output_file)
							prog_data[ "recoms" ].append( event )

							bir += 1

						# Other recombinations have come before, meaning this IS a recombination event
						# There's just no recom_bool to trigger it this time
						else: 
							recom_count += 1
							prog_recom_count += 1
							if arm_side == "left":
								event_size = abs( recom_frag_dict[ "positions"][0] )
							else: 
								event_size = abs( recom_frag_dict[ "positions"][0] - chr_lengths[ chr_name ])
							event = ",".join( [chr_name, str(recom_frag_dict[ "positions"][0]), str(event_size)] )
							print( prog_name + "," + event, file = output_file)
							prog_data[ "recoms" ].append( event )

							bir += 1

				else: 
					if first_recom:
						if recom_bool: 

							# First no longer applies, so turn it off! 
							first_recom = 0

							# recom_frag_dict is reinitiated down here in the normal cases
							recom_frag_dict[ "positions" ] = [ posit_list[ index ] ]
							recom_frag_dict[ "markers" ] = [ marker_list[ index ] ]

							# Turn these on if you want to keep the full potential length of the first recombination event
							# recom_frag_dict[ "positions" ] = [ posit_list[ stor_idx ] ]
							# recom_frag_dict[ "markers" ] = [ marker_list[ stor_idx ] ]

							# Reassign stor_idx for the 
							stor_idx = index

					# Average, everyday recombination events: 
					else: 
						if recom_bool: 
							
							# Do the regular recombination stuff again
							recom_count += 1
							prog_recom_count += 1
							event_size = abs(recom_frag_dict[ "positions"][0] - recom_frag_dict[ "positions"][-1])
							event = ",".join( [chr_name, str(recom_frag_dict[ "positions"][0]), str(event_size)] )
							print( prog_name + "," + event, file = output_file)
							prog_data[ "recoms" ].append( event )

							# Reinitiate the recom_frag_dict at the current position
							recom_frag_dict[ "positions" ] = [ posit_list[ index ] ]
							recom_frag_dict[ "markers" ] = [ marker_list[ index ] ]
							stor_idx = index

		print( prog_name + "," + str(prog_recom_count), file = count_file)
	
	print( "Total progeny: ", total_progeny )
	print( "Total recombinations: ", recom_count )
	print( "Number of putative break-induced replications: ", bir )

	return recombs, total_progeny # Recombs is just a dictionary with progeny names as keys and then a list of strings with recombination "events" 
	# Sample event: event = ",".join( [chromosomes[i], str( left_side_centro ), str(event_size)] )

	f2_file.close()
	output_file.close()
	count_file.close()


def recom_Stacker( chr_lengths, recom_dict, total_progeny, output_file_name, resolution = 100000, \
					a = "1", b = "2", h = "n", missing = "-", strange = "x" ): 
	
	"""Creates a dictionary of chromosomes and bin sizes that wll be used for histograms to depict recombination hotspot for each allele type."""   

	# Compute the universal position for each chromosome into the total variable. 
	# This is to translate marker positions into universal positions for graphing 

	total = 0
	univ_chr_dict = {}
	for chromosomes, lengths in chr_lengths.items(): 
		univ_chr_dict[ chromosomes ] = total
		total += int(lengths+resolution)

	# Build a dictionary to store the recombination events and for tallying the markers
	recom_bins = {} 
	tally_dict = {} 

	for chromo, length in chr_lengths.items(): 
		chr_bin_dict = {}
		tal_bin_dict = {}

		stop = length-(length%resolution)+resolution
		
		for i in range(0, stop, resolution):
			chr_bin_dict[ i ] = 0
			tal_bin_dict[ i ] = { a: 0, b: 0, h: 0, missing: 0, strange: 0 }

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
			start_bin = position - (position % resolution)  
			
			# Put that into the appropriate bin of tally_dict:
			recom_bins[ chromosome ][ start_bin ] += 1

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
	header = "chr,chr.pos,univ.pos,a,b,h,s,missing,recoms,recom.by.col,recom.by.mkr"
	  
	print( header, file = recom_out_file )
	for chromosome, chr_bins in recom_bins.items(): 
		for bins, recom_freq in chr_bins.items(): 
			# I'm gonna call the tallies using the recombination data because the positions are identical and I want them in the same file. 
			prt_pos = str(bins)
			prt_univ_pos = str(bins + univ_chr_dict[ chromosome ])

			int_a = tally_dict[ chromosome ][ bins ][ a ]
			int_b = tally_dict[ chromosome ][ bins ][ b ]
			int_h = tally_dict[ chromosome ][ bins ][ h ]
			int_s = tally_dict[ chromosome ][ bins ][ strange ]
			int_m = tally_dict[ chromosome ][ bins ][ missing ]
			int_freq = recom_freq

			total = int_a + int_b + int_h + int_s + int_m # Every marker, present or absent

			mkr_loci_count = total // total_progeny
			id_mkr_cnt = total - int_m # Every marker with an identifiable genotype (excludes missing markers)

			try: 
				recom_by_loci_cnt = round( int_freq / mkr_loci_count, 4 )
				recom_by_mkr_cnt = round( int_freq / id_mkr_cnt, 4 )

			except ZeroDivisionError: 
				recom_by_loci_cnt = 0
				recom_by_mkr_cnt = 0 

			prt_lc = str( recom_by_loci_cnt )
			prt_rmc = str( recom_by_mkr_cnt )

			prt_a, prt_b, prt_h, prt_s, prt_m, prt_freq = str(int_a), str(int_b), str(int_h), str(int_s), str(int_m), str(int_freq)

			print( ",".join( [ chromosome, prt_pos, prt_univ_pos, prt_a, prt_b, prt_h, prt_s, prt_m, prt_freq, prt_lc, prt_rmc ] ), file = recom_out_file )
			
	recom_out_file.close()

	return recom_bins, resolution


def centromerely_math( chr_lengths, centromere_locations, recom_dict, output_file_name, \
						resolution = 100000, count_incentro_bool = 0, left_arm = "left",\
						right_arm = "right", centromere = "centromere" ): 

	""" Calculate cumulative recombination as a function of distance from centromeres """

	# Compute the universal position for each chromosome into the total variable. 
	# This is to translate marker positions into universal positions for graphing.

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
			# recom_side = recom_start + recom_size
			recom_side = recom_start
			
			if chro not in centro_recoms:  # Using recom_side because that's what I use to count recombination for each position
				centro_recoms[ chro ] = { left_arm:{}, right_arm:{}, centromere:{} }
				sorted_recoms[ chro ] = { left_arm:[], right_arm:[], centromere:[] } # This dictionary and the next are for later
				cyuum_recoms[ chro ]  = { left_arm:{}, right_arm:{}, centromere: 0 }


			# Make sure the recombination isn't anywhere else
			put_er_in_bool = (recom_side not in centro_recoms[ chro ][ left_arm   ]) and \
			                 (recom_side not in centro_recoms[ chro ][ right_arm  ]) and \
			                 (recom_side not in centro_recoms[ chro ][ centromere ])

			# If the recombination can be put into the dictionary, then put it in there. 
			if put_er_in_bool: 
				if recom_side <= centromere_locations[ chro ][0]: 
					centro_recoms[ chro ][right_arm][ recom_side ] = 1
				elif recom_side >= centromere_locations[ chro ][1]: 
					centro_recoms[ chro ][left_arm][ recom_side ] = 1
				else: 
					centro_recoms[ chro ][ centromere ][ recom_side ] = 1

			# If it's already in there, then add another event to it
			else:
				if recom_side <= centromere_locations[ chro ][0]: 
					centro_recoms[ chro ][right_arm][ recom_side ] += 1
				elif recom_side >= centromere_locations[ chro ][1]: 
					centro_recoms[ chro ][left_arm][ recom_side ] += 1
				else: 
					centro_recoms[ chro ][ centromere ][ recom_side ] += 1

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


def main(create_tester = 0, chr_file_name = "calbicans_chromosomes.csv"):

	"""Main combines all the functions into one "smooth" workflow!"""

	arg_list = []
	for arg in sys.argv:
		arg_list.append( arg )
	
	path_and_name = arg_list[1] # Parse the path and file name
	# use_Strangers = int(arg_list[2])

	# Split the name into two pieces: the directory path and the file name itself
	dir_path, file_name = parse_path( path_and_name )

	# Obtaining the chromosomal information from the file!
	try: 
		# chr_file = open( dir_path + chr_file_name, "r" ) # Chromosome files in subfolders; you may need to change chr_file_name argument to use the file you want
		chr_file = open( chr_file_name, "r" )
	except IOError: 
		print( "Couldn't find the chromosomes file:\n " + chr_file_name + ". \nQuitting..." )

	centromere_locations = {}
	chr_lengths = {}
	chr_labs = {}
	header_skip = 1
	for lines in chr_file: 
		if header_skip: 
			header_skip = 0 
			continue
		line_list = lines.split(",")
		chr_lengths[ line_list[0] ] = int(line_list[1])
		centromere_locations[ line_list[0] ] = [int(line_list[2]),int(line_list[3])]
		chr_labs[line_list[0]] = line_list[4].strip()

	# If performing a test run, this will be the section that is activated
	test_bool = path_and_name.strip().lower() == "test"
	if test_bool: 
		path_and_name = "test/test.csv"
		dir_path, file_name = parse_path( path_and_name )
		test_creator( path_and_name, 65, 1200, 4)

	marker_cleaner( dir_path, file_name )
	
	cleaned_file_name = dir_path + "cleaned_" + file_name
	undesirable_file_name = dir_path + "undesirables_" + file_name

	four_to_F2( cleaned_file_name )

	f2_file_name = "f2_cleaned_" + file_name
	ploid_f2_file_name = "f2_recleaned_" + file_name
	if test_bool:
		blacklist_file_name = "test/blk_list_test.csv"
	else:
		# The blacklist file will only remove items that are present, so I've added anything I want to remove from either cross
		blacklist_file_name = "progeny_ploidy_blacklist.csv"
		# blacklist_file_name = "blank_blacklist.csv" # Use this line to compare "before and after" blacklisting
	
	deploid_dict, number_o_progeny = deploidy( dir_path, f2_file_name, ploid_f2_file_name, chr_labs, blacklist_file_name )
	stack_bar_file_name = dir_path + "stackbar_" + file_name 
	marker_tally_ho( deploid_dict, stack_bar_file_name )
	os.system( "Rscript recombination_bins.R " + dir_path + chr_file_name + " " + stack_bar_file_name + " s " )

	prog_recom_count_filename = dir_path + "recom-count_" + file_name
	stats_file_name = dir_path + "statistics_" + file_name
	recom_dict, total_progeny = recombination_analysis( dir_path, ploid_f2_file_name, stats_file_name, prog_recom_count_filename, centromere_locations, chr_lengths )
	os.system( "Rscript recombination_bins.R " + dir_path + chr_file_name + " " + stats_file_name + " h " + prog_recom_count_filename  )

	resolve_this = 10000
	recom_data_name = dir_path + "recom_" + file_name
	recom_bins, resolution = recom_Stacker( chr_lengths, recom_dict, total_progeny, recom_data_name, resolution = resolve_this ) # Stack recombinations
	os.system( "Rscript recombination_bins.R " + dir_path + chr_file_name + " " + recom_data_name + " g" )

	centro_data_name = dir_path + "centro_" + file_name
	centromerely_math( chr_lengths, centromere_locations, recom_dict, centro_data_name, resolution = resolve_this ) # takes The dictionary from recombination_analysis(), not recom_Stacker()
	os.system( "Rscript recombination_bins.R " + dir_path + chr_file_name + " " + centro_data_name + " c" )

main()