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


def marker_cleaner(sample_list, file_name): 
	"""This function removes markers that are missing from eiter parent and markers that are
	ambiguous between parents (the parents share an allele). """

	keepers = []
	chuckers = []

	length = len(sample_list[0])

	chromosomes = sample_list[0]
	chr_pos     = sample_list[1]
	parent1     = sample_list[2] #SC5314
	parent2     = sample_list[3] #529L or P60002

	for index in range(1,length): # Skip the header; start at 1

		try: 
			p1b = parent1[index][0] + parent1[index][-1]
			p2b = parent2[index][0] + parent2[index][-1]
		except IndexError: 
			print( "Something is wrong with the parental genotypes: ", chromosomes[index], ": ", chr_pos[index] )
			print( "Parent 1 (SC): ", p1b )
			print( "Parent 2 (529L/P6): ", p2b )

		try:
			chr_pos[index] = int( chr_pos[index] ) # Convert to ints so they can get sorted properly
		except IndexError: 
			print( "Something is wrong with the chromosome position: ", chromosomes[index], ": ", chr_pos[index] )
			quit()

		# if they share an allele! Make sure I know...
		chucker_bool = (p1b[0] == p2b[0]) or (p1b[0] == p2b[1]) or (p1b[1] == p2b[0]) or (p1b[1] == p2b[1])

		if chucker_bool:
			chuckers.append( index )
		else: 
			keepers.append( index )

	# Convert the data to 4way to just put it through the pipeline! It'll be much simpler that way.  


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


def main():

	arg_list = []
	for arg in sys.argv:
		arg_list.append( arg )
	
	path_and_name = arg_list[1] 

	dir_path, file_name = parse_path( path_and_name )

	data_file = open( file_name, "r" )

	trp_data_list = csv_reader( data_file )

	data_list = transpose( trp_data_list )

	chromosome    = data_list[0]
	chr_pos       = data_list[1] 
	mutation_type = data_list[2]
	REF           = data_list[3] # The reference genotype (SC5314)
	ALT           = data_list[4] # The parental variant phenotype (either 529L or P60002, depending on the strain)
	hom_ref       = data_list[5] # Number of samples that are homozygous for the reference
	hom_alt       = data_list[6] # Number of samples that are homozygous for the variant
	het_genos     = data_list[7] # Number of samples that have heterozygous genotypes
	het_var       = data_list[8] # Number of samples with 2 alleles that don't belong to SC5314 - both are from the other parent
	
	# Genotypes for each parent and progeny
	SC5314                 = data_list[9]
	five29L                = data_list[10]
	P60002                 = data_list[11]
	MAY_103_529L_tet       = data_list[12] # SCx529L Tetraploid
	MAY_316_P600_tet       = data_list[13] # SCxP6 Tetraploid
	MAY_154_529L_prog      = data_list[14] # SCx529L Progeny
	MAY_155_529L_prog      = data_list[15] # SCx529L Progeny
	MAY_156_529L_prog      = data_list[16] # SCx529L Progeny
	MAY_157_529L_prog      = data_list[17] # SCx529L Progeny
	MAY_158_529L_prog      = data_list[18] # SCx529L Progeny
	MAY_332_P600_prog      = data_list[19] # SCxP60002 Progeny
	MAY_333_P600_prog      = data_list[20] # SCxP60002 Progeny

	SCx529L = [ chromosome, chr_pos, SC5314, five29L, MAY_103_529L_tet, MAY_154_529L_prog, MAY_155_529L_prog, MAY_156_529L_prog, MAY_157_529L_prog, MAY_158_529L_prog ]
	SCxP6   = [ chromosome, chr_pos, SC5314, P60002, MAY_316_P600_tet, MAY_332_P600_prog, MAY_333_P600_prog ]

	print( "SC5314 x 529L" )
	marker_cleaner( SCx529L, "cleaned_SCx529L.csv" )
	print() 
	print( "SC5314 x P60002" )
	marker_cleaner( SCxP6, "cleaned_SCxP6.csv" )

main()