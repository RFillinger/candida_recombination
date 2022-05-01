

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
				max_dist = 100, start_dist = 1, increment = 1, total_tests = 10000 ):

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

		srtble_recom = [ progeny_name, chromosome, break_pos, break_dist ]
		empty_list.append( srtble_recom )

	sorted_list = sorted( empty_list, key = operator.itemgetter(1,2))

	coin_flip_list = list(range(1,101))
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

			coin_flip = random.choice( coin_flip_list )
			if coin_flip >= 51: 
				towards_centromere = 1
			else: 
				towards_centromere = 0

			# Find where this is relative to centromere
			rando_event = random.choice( sorted_list )
			index = sorted_list.index(rando_event)

			progeny_name = rando_event[0]
			chromosome   = rando_event[1]
			break_pos    = rando_event[2] 
			break_dist   = rando_event[3]

			if break_pos >= centromere_locations[ chromosome ][1]: 
				
				arm = "right"

				if towards_centromere: 
				
					try: 
						previous_event = sorted_list[index-1]

						# Check to see if it crosses the centromere. 
						if (break_pos - distance) < centromere_locations[chromosome][1]: 
							continue

						if chromosome != previous_event[1]: 
							continue

						# Check to see if there is a successfully found recombination
						success_bool = abs(break_pos-previous_event[2]) <= distance
						if success_bool: 
							successful_forays += 1
						test_number += 1

					except IndexError: 
						pass # There is no recombination found; this is the end of the chromosome No increment of test_number
				

				else: 
					if break_dist > distance: 
						test_number += 1 # Nothing happens, no recom was found with this trial
						
					else: 
						successful_forays += 1
						test_number += 1


			elif break_pos <= centromere_locations[ chromosome ][0]: 
				
				arm = "left"

				if not towards_centromere: # This was backwards and doesn't work without the "not"... sorry!
					if break_dist > distance: 
						test_number += 1 # Nothing happens, no recom was found with this trial
						
					else: 
						successful_forays += 1
						test_number += 1
				else:
					try: 
						previous_event = sorted_list[index+1]

						# Check to see if it crosses the centromere. 
						if (break_pos + distance) > centromere_locations[chromosome][0]: 
							continue

						if chromosome != previous_event[1]: 
							continue

						# Check to see if there is a successfully found recombination
						success_bool = abs(break_pos-previous_event[2]) <= distance
						if success_bool: 
							successful_forays += 1
						test_number += 1

					except IndexError: 
						pass # There is no recombination found; this is the end of the chromosome No increment of test_number

			else: 
				arm = "centromere"
				# print(rando	

		probability_dict[distance] = (round((successful_forays/(test_number*len(unique_progeny))),5))		

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