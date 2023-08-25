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
			try: 
				row.append(item[i])
			except IndexError: 
				print( item, i  )
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


def snp2mkr( SNP, cat_ID ):

	chromosome = SNP[0]
	
	try: 
		chr_pos = SNP[1]
	except IndexError: 
		print( SNP )

	mkr_list = []

	p1_raw = SNP[0][0] + SNP[0][-1]
	p2_raw = SNP[1][0] + SNP[1][-1]

	# Convert the data to 4way here
	geno_list = ["a", "b", "c", "d"]

	empty_set = set()
	for genos in SNP:

		empty_set.add(genos[0])
		empty_set.add(genos[-1])


	if len(empty_set) > 4:
		print(SNP) 
		print( empty_set )
		print( "There are 5 alleles at the marker:", chromosome, chr_pos)
		print( "This pipeline can't really deal with that... so it's quitting. ")
		quit()

	marker_dict = {}
	for idx, elements in enumerate(empty_set):
		marker_dict[ elements ] = geno_list[idx]

	new_mkr_list = []
	for mkrs in SNP:

		new_mkr = ""
		for char in mkrs: 

			if (char == "/") or (char == "|"): 
				continue

			if (char == ".") or (char == "*"): 
				new_mkr = "-"
				break
			
			else:
				new_mkr = new_mkr + marker_dict[char]

		new_mkr_list.append( new_mkr )

	return new_mkr_list  


def this_is_the_way():

	"""This function converts "diploid" SNP data FROM 529L/P60002 parents and progeny 
		into 4way data for use in recombination_analysis.py """

	arg_list = []
	for arg in sys.argv:
		arg_list.append( arg )
	
	path_and_name = arg_list[1] 
	dir_path, file_name = parse_path( path_and_name )

	if "test" in file_name.lower(): 
		test = 1
	else: 
		test = 0
		
	data_file = open( file_name, "r" )
	data_list = csv_reader( data_file )

	five29L_cross = []
	p60002_cross = []

	header = 1
	for index, snps in enumerate(data_list): 

		# print( snps )

		chromosome        = snps[0]
		chr_pos           = snps[1]
		par_1            = snps[9]
		five29L           = snps[10]
		P60002            = snps[11]
		MAY_103_529L_tet  = snps[12] # SCx529L Tetraploid
		MAY_316_P600_tet  = snps[13] # SCxP6 Tetraploid
		MAY_154_529L_prog = snps[14] # SCx529L Progeny
		MAY_155_529L_prog = snps[15] # SCx529L Progeny
		MAY_156_529L_prog = snps[16] # SCx529L Progeny
		MAY_157_529L_prog = snps[17] # SCx529L Progeny
		MAY_158_529L_prog = snps[18] # SCx529L Progeny
		MAY_332_P600_prog = snps[19] # SCxP60002 Progeny
		MAY_333_P600_prog = snps[20] # SCxP60002 Progeny

		cat_ID = str(index)

		if "Ca19-mtDNA" in chromosome: 
			continue

		if header: 
			header = 0
			# five29L_cross.append( ["Chromosome", "Chr_Position", "# Catalog ID", snps[9], snps[10], snps[12]] + snps[14:19] )
			five29L_cross.append( ["Chromosome", "Chr_Position", "# Catalog ID", snps[9], snps[10]] + snps[14:19] )
			p60002_cross.append( ["Chromosome", "Chr_Position", "# Catalog ID", snps[9], snps[11], snps[13]] + snps[19:] )
			continue

		# SCx529L = [ par_1, five29L, MAY_103_529L_tet, MAY_154_529L_prog, MAY_155_529L_prog, \
		# 											MAY_156_529L_prog, MAY_157_529L_prog, MAY_158_529L_prog ] # Include the tetraploid

		SCx529L = [ par_1, five29L, MAY_154_529L_prog, MAY_155_529L_prog, MAY_156_529L_prog, MAY_157_529L_prog, MAY_158_529L_prog ] # Exclude tetraploid

		SCxP6   = [ par_1, P60002, MAY_316_P600_tet, MAY_332_P600_prog, MAY_333_P600_prog ]

		lists = [SCx529L, SCxP6]

		new_mkr = []
		for cross in lists: 
			new_mkr.append(snp2mkr(cross, cat_ID))

		new_529L_mkr = [chromosome, chr_pos, cat_ID] + new_mkr[0]
		new_P6_mkr = [chromosome, chr_pos, cat_ID] + new_mkr[1]

		five29L_cross.append( new_529L_mkr )
		p60002_cross.append( new_P6_mkr )

	if test: 

		new_5_file = open( "529L_TEST_4way.csv", "w" )
		# csv_printer( transpose(five29L_cross), new_5_file )
		csv_printer( five29L_cross, new_5_file )

		new_P6_file = open( "P600_TEST_4way.csv", "w" )
		# csv_printer( transpose(p60002_cross), new_P6_file )
		csv_printer( p60002_cross, new_P6_file )


	else: 

		new_5_file = open( "wild_isolates_recombination.csv", "w" )
		# csv_printer( transpose(five29L_cross), new_5_file )
		csv_printer( five29L_cross, new_5_file )

		# new_P6_file = open( "ltd_P6_4way.csv", "w" )
		# csv_printer( transpose(p60002_cross), new_P6_file )
		# csv_printer( p60002_cross, new_P6_file )


def main():

	"""This function converts "diploid" SNP data FROM 529L/P60002 parents and progeny 
		into 4way data for use in recombination_analysis.py """

	arg_list = []
	for arg in sys.argv:
		arg_list.append( arg )
	
	path_and_name = arg_list[1] 
	dir_path, file_name = parse_path( path_and_name )

	data_file = open( path_and_name, "r" )
	data_list = csv_reader( data_file )

	five29L_cross = []

	header = 1
	for index, snps in enumerate(data_list): 

		chromosome = snps[0]
		chr_pos    = snps[1]
		par_1      = snps[9]
		par_2      = snps[10]
		prog_1     = snps[11]
		prog_2 	   = snps[12]

		if "1" in par_1: 
			print( snps )

		cat_ID = str(index)

		if "Ca19-mtDNA" in chromosome: 
			continue

		if header: 
			header = 0
			five29L_cross.append( ["Chromosome", "Chr_Position", "# Catalog ID", snps[9], snps[10], snps[11], snps[12]] )
			continue

		cross = [ par_1, par_2, prog_1, prog_2 ]

		new_mkr = snp2mkr(cross, cat_ID)

		new_529L_mkr = [chromosome, chr_pos, cat_ID] + new_mkr

		five29L_cross.append( new_529L_mkr )

	new_5_file = open( dir_path + "WIR2.csv", "w" )

	csv_printer( transpose(five29L_cross), new_5_file )

main()
