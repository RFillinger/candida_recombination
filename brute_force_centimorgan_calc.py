

import operator
import random
import sys


def parse_path( path_n_name, sep = "/" ): 

	''' Takes a file path and name and splits it into two variables to be used for output. '''
	path_list = path_n_name.split( sep )
	file_name = path_list[-1]
	dir_path = "/".join(path_list[:-1]) + "/"
	if dir_path == "/": 
		dir_path = ""

	return dir_path, file_name


def csv_reader( input_file ):

	'''Reads in .csv files into a list'''

	line_list = []
	for lines in input_file: 
		line_list.append( lines.strip().split(",") )

	return line_list 


def univ_pos( chr_lengths ):

	total = 0 
	universal_positions = {}
	for names, length in chr_lengths.items(): 
		total = total + length 
		universal_positions[ names ] = total
	
	return universal_positions


def centimomos( dir_path, file_name, chr_file_name = "calbicans_chromosomes.csv", 
				max_dist = 850, start_dist = 750, increment = 1, total_tests = 100000 ):

	"""Brute forcing centimorgan's calculation."""

	try: 
		# chr_file = open( dir_path + chr_file_name, "r" ) # Chromosome files in subfolders; you may need to change chr_file_name argument to use the file you want
		chr_file = open( dir_path + chr_file_name, "r" )
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

	recom_file = open( dir_path + file_name, "r" )

	empty_list = []
	unique_progeny = set()
	chromosome_set = set()
	
	header = 1 
	for recom in recom_file: 
		if header: 
			header = 0 
			continue

		recom_list = recom.split(",")

		progeny_name = recom_list[0]
		chromosome   = recom_list[1]
		break_pos    = int(recom_list[2]) # Position of a detected break
		break_dist   = int(recom_list[3]) # Distance to the next detectable break (distance from break to next break FURTHER AWAY FROM CENTROMERE)

		unique_progeny.add(progeny_name) # Number of unique progeny for normalizing recombination probability over these events
		chromosome_set.add(chromosome)

		srtble_recom = [ progeny_name, chromosome, break_pos, break_dist ]
		empty_list.append( srtble_recom )

	chromosome_list = list(chromosome_set)

	sorted_list = sorted( empty_list, key = operator.itemgetter(1,2))

	
	recom_chr_dict = {}
	for recoms in sorted_list: 

		progeny_name = recoms[0]
		chromosome   = recoms[1]
		break_pos    = recoms[2] # Position of a detected break
		break_dist   = recoms[3]

		new_event = [progeny_name, break_pos, break_dist]
		
		if chromosome in recom_chr_dict: 
			recom_chr_dict[chromosome].append( new_event )
		else: 
			recom_chr_dict[chromosome] = []

	probability_dict = {}

	print("\nProgress: ")
	distance = start_dist
	while distance <= max_dist:

		perc_trip = round( 100*(distance - start_dist)/(max_dist - start_dist), 0 )

		if perc_trip % increment == 0: 
			print( str(int(perc_trip)) + "%" )

		successful_forays = 0

		test_number = 1
		while test_number <= total_tests: 

			# Do all the shit in here. 
			rando_chr = random.sample(chromosome_list, 1)[0]
			rando_start = random.randint(1,chr_lengths[rando_chr])
			rando_end = rando_start + distance

			if rando_end > chr_lengths[rando_chr]:
				rando_start = chr_lengths[rando_chr] - distance 
				rando_end = chr_lengths[rando_chr]

			for recoms in recom_chr_dict[ rando_chr ]: 
				break_pos = recoms[1]
				
				if (rando_start <= break_pos) and (rando_end >= break_pos): 
					successful_forays += 1
					break

			test_number += 1

		probability_dict[distance] = round((successful_forays/(total_tests*len(unique_progeny))),5)

		distance += increment

	print( "Done!\n")

	return probability_dict


def main():

	"""Main combines all the functions into one "smooth" workflow!"""

	arg_list = []
	for arg in sys.argv:
		arg_list.append( arg )
	
	# Parse the path and file name
	path_and_name = arg_list[1] 

	# Split the name into two pieces: the directory path and the file name itself
	dir_path, file_name = parse_path( path_and_name )

	probabilities = centimomos( dir_path, file_name )

	print(probabilities)

main()