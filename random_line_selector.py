import random
import sys


def parse_path( path_n_name, sep = "/" ): 

	''' Takes a file path and name and splits it into two variables to be used for output. '''
	path_list = path_n_name.split( sep )
	file_name = path_list[-1]
	dir_path = sep.join(path_list[:-1]) + sep
	if dir_path == sep: 
		dir_path = ""

	return dir_path, file_name


def csv_reader( input_file, sep = "," ):

	'''Reads in .csv files into a list'''

	line_list = []
	for lines in input_file: 
		line_list.append( lines.strip().split(sep) )

	return line_list


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


def main(): 

	arg_list = []
	for arg in sys.argv:
		arg_list.append( arg )
	
	# Command should read something like: python3 random_line_selector.py path/to/file/file.txt 10
	path_and_name = arg_list[1]
	try: 
		number_o_lines = int(arg_list[2])
	except TypeError: 
		print("Second argument needs to be a number less than the length of the file. Quitting.")

	path, file_name = parse_path(path_and_name)

	target_file = open( path_and_name, "r" )

	file_list = csv_reader( target_file )

	rando_line_indices = set()
	while len(rando_line_indices) < number_o_lines: 
		rando = random.randint(0,len(file_list)-1)
		if rando == 0: # Avoid the header
			continue
		rando_line_indices.add(rando)

	rando_list = list(rando_line_indices)
	
	new_line_list = []
	for indices in rando_list: 
		new_line_list.append(file_list[indices])

	output_file = open( path + "rando_recombs_" + file_name, "w" )
	csv_printer(new_line_list, output_file)

main()