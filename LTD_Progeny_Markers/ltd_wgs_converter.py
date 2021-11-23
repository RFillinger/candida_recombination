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
	SC5314        = data_list[9]
	five29L       = data_list[10]
	P60002        = data_list[11]
	MAY_103       = data_list[12]
	MAY_316       = data_list[13]
	MAY_154       = data_list[14]
	MAY_155       = data_list[15]
	MAY_156       = data_list[16]
	MAY_157       = data_list[17]
	MAY_158       = data_list[18]
	MAY_332	      = data_list[19]
	MAY_333       = data_list[20]

	# for lines in data_list: 
	# 	print( lines[0] )

	# Next thing to do is to convert the parental genotypes into 1's, 2's, and n's

main()