
import sys
import numpy
import statistics
import operator
import random
import copy
import numpy

# Finds a couple things: 
# 1) Incorporate the location of recombination events into the data

def csv_reader( input_file ):

	'''Reads in .csv files into a list'''

	line_list = []
	for lines in input_file: 
		line_list.append( lines.strip().split(",") )

	return line_list 


def transpose(l1):

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


def element_swap( listus, item_is, item_should_be ): 

	for n, i in enumerate(listus): # Change chromosome 8 to "R"
		if i == item_is: 
			listus[n] == item_should_be

	return listus


def csv_printer( csv_list, output_file ):


	for items in csv_list: 

		line = ",".join(items)

		print( line, file = output_file)


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


def test_creator(): 

	shoot_trouble = 0
	number_of_prog_to_generate = 66
	length = 200

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

			if num == 1: # only "aa" is available
				
				leel_list.append( prog_combs[0] ) 

			elif ( num > 1 ) and ( num <= 6 ): # Only "a" or "b" allele combinations are possible

				rint_6 = random.randint(1, 3)
				leel_list.append( prog_combs[rint_6-1])

			elif ( num > 6 ) and ( num <= 10 ): # A, B, and C alleles are available

				rint_10 = random.randint(1, 6)
				leel_list.append( prog_combs[rint_10-1])

			else: # All aleles A through D are available

				rint_11 = random.randint(1, 10)
				leel_list.append( prog_combs[rint_11-1])

		big_boi.append( leel_list )

		index += 1

	if shoot_trouble: 
		print( big_boi )
		print( dist_list )
		print( type( big_boi) )

	header_parts = albicans_header_creator(length)

	type_swap_header_list = [ [], [], [] ]
	
	for n, sublists in enumerate( header_parts ): 

		for elements in sublists: 

			type_swap_header_list[n].append( str(elements) )


	test_file = open( "recombination_test.csv", "w" )

	print( ",".join(type_swap_header_list[0]), file = test_file )
	print( ",".join(type_swap_header_list[1]), file = test_file )
	print( ",".join(type_swap_header_list[2]), file = test_file )


	csv_printer( big_boi, test_file, transpose_bool = 1 ) # Change this to print the file and we're good to go. 
	test_file.close()


def marker_cleaner(): 

	"""This function removes markers that are missing from eiter parent and markers that are
	ambiguous between parents (the parents share an allele). """

	print_undesirables = 1

	arg_list = []

	for arg in sys.argv:
	    arg_list.append( arg )

	four_way_file = open( arg_list[1], "r" )

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
	new_file = open( "cleaned_" + arg_list[1], "w" )
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
		new_un_file = open( "undesirables_" + arg_list[1], "w" )
		for lines in trp_chuckers: 
			print( ",".join(lines), file = new_un_file )
		new_un_file.close()


def four_to_F2(): 

	"""Converts 4-way alele data to F2 
	ALL INPUT INTO THIS FUNCTION MUST GO THROUGH THE 'marker_cleaner' FUNCTION FIRST!!!"""

	arg_list = []

	for arg in sys.argv:
		arg_list.append( arg )

	cleaned_marker_file = open( arg_list[1], "r" )

	untrp_cleaned_markers = csv_reader( cleaned_marker_file )

	cleaned_markers = transpose(untrp_cleaned_markers)

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

	new_file = open( "f2_" + arg_list[1], "w")
	csv_printer( trp_f2, new_file )
	new_file.close()


def recombination_analysis():

	arg_list = []

	for arg in sys.argv:
		arg_list.append( arg )

	f2_file_name = arg_list[1]
	output_name  = arg_list[2] 

	f2_file = open( f2_file_name, "r" )
	big_list = csv_reader( f2_file )

	output_file = open( output_name, "w" )
	print( "Sample Name, Filled Markers, Parent 1 Shared Markers, Parent 2 Shared Markers, Neither Parent Markers, Missing Markers, Number of Recombination Events, Estimated Recombination Site Length (bp), St. Dev. of Recombination Site Length", file = output_file)

	mrkr_list = big_list[0] #big_list is a 2D array of the .csv file 
	chromosome = big_list[1]
	chr_pos = big_list[2]

	parent_1 = big_list[3]
	parent_2 = big_list[4]

	prog_allele_list = big_list[5:] #This is all the progeny stored in the 2D array

	m1 = "1"
	m2 = "2"
	m3 = "n"

	for progeny in prog_allele_list:

		prog_name = progeny[0]

		recomb_lengths = []

		temp_genotype = progeny[1]
		temp_position = chr_pos[1]
		temp_chr = chromosome[1]

		mean_dis = "-" # These get turned into ints later, the "-" is for data labelling that occurs during an error
		stdev_dis = "-"

		num_1 = 0
		num_2 = 0
		num_n = 0 
		total = 0
		missing = 0

		if progeny[0] == m1:
			num_1 += 1 
			total += 1 
		elif progeny[0] == m2:
			num_2 += 1
			total += 1 
		elif progeny[0] == m3:
			num_n += 1 
			total += 1
		else:
			missing += 1

		index = 2 # Start in the data; not the sample name

		while index < len( mrkr_list ):

			if progeny[ index ] == m1: # Gathering data for the output
				num_1 += 1 
				total += 1 
			elif progeny[ index ] == m2:
				num_2 += 1
				total += 1 
			elif progeny[ index ] == m3:
				num_n += 1 
				total += 1 
			else:
				missing += 1

			if chromosome[ index ] != temp_chr: # If the chromosome changes
				
				if ( progeny[ index ] != "-" ):

					temp_chr = chromosome[ index ] 
					temp_position = chr_pos[ index ]

					temp_genotype = progeny[ index ]

					index += 1
					continue

				else: 
					index += 1
					continue

			good_to_go = 0 # This is a little convoluted but it works; it's looking for an exception when parental alleles are identical but can be differentiated in progeny

			if ( parent_1[ index ] == parent_2[ index ] ) and \
			   ( progeny[ index ] != parent_1[ index ] ):

			   good_to_go = 1

			if ( progeny[ index ] != temp_genotype )      and \
			   ( progeny[ index ] != "-" )                and \
			   (( parent_1[ index ] != parent_2[ index ] ) or \
			   good_to_go):

			   # If progeny genotype changed from left to right AND
			   # Progeny genotype is not empty                  AND
			   # ( Parent_1's allele set is different from Parent_2's OR 
			   # Parent_1 and Parent_2's genotype are the same BUT Progeny genotype is different from either of them (that's the good_to_go boolean))
				
				try:
					segment_len = int( chr_pos[ index ] ) - int( temp_position )	

					recomb_lengths.append( segment_len )

					temp_position = chr_pos[ index ]
					temp_genotype = progeny[ index ]
				
				except ValueError: 
					pass
			
			index += 1

		try: # Too few data points to actually get any information; error will trip if there is nothing in the list

			mean_dis = int( round( statistics.mean( recomb_lengths ), 0))
			stdev_dis = int( round( statistics.stdev( recomb_lengths ), 0))

			print( prog_name + "," + str( total ) + "," + str( num_1 ) + "," + str( num_2 ) + "," + str( num_n ) + "," + str( missing ) + "," + \
			       str( len( recomb_lengths )) + "," + str( mean_dis ) + "," + str( stdev_dis ), file = output_file )

		except statistics.StatisticsError: 

			print( prog_name + "," + str( total ) + "," + str( num_1 ) + "," + str( num_2 ) + "," + str( num_n ) + "," + str( missing ) + "," + \
				   str( len( recomb_lengths )) + "," + "-" + "," + "-", file = output_file )

	f2_file.close()
	output_file.close()


# test_creator()
# marker_cleaner()
# four_to_F2()
recombination_analysis()


