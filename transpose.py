import sys

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


def parse_path( path_n_name, sep = "/" ): 

	''' Takes a file path and name and splits it into two variables to be used for output. '''
	path_list = path_n_name.split( sep )
	file_name = path_list[-1]
	dir_path = "/".join(path_list[:-1]) + "/"
	if dir_path == "/": 
		dir_path = ""

	return dir_path, file_name


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
	
	path_and_name = arg_list[1]

	dir_path, file_name = parse_path( path_and_name )
	file = open( path_and_name, "r" )
	array = csv_reader( file )
	trp_array = transpose( array )
	new_file = open( dir_path + "trp_" + file_name, "w" )
	csv_printer(trp_array, new_file )


main()