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


def csv_printer( csv_list, output_file, sep = ",", include_keys = 0 ):

	if isinstance(csv_list, list): 
		for items in csv_list: 
			line = sep.join(items)
			print( line, file = output_file)
	elif isinstance(csv_list, dict):
		if include_keys:
			for keys, values in csv_list.items(): 
				line = sep.join([keys]+values)
				print( line, file = output_file)
		else: 
			for keys, values in csv_list.items(): 
				line = sep.join(values)
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


def snp2mkr( SNP, cat_ID, ref_allele, alt_allele, pass_on_5_mkrs = True ):

	new_snp = []
	for genotypes in SNP: 

		new_geno = genotypes.replace(".", ref_allele, genotypes.count("."))

		new_snp.append(new_geno)

	SNP = new_snp

	mkr_list = []
	
	# Convert the data to 4way here
	geno_list = ["a", "b", "c", "d"]

	empty_set = set()
	for genos in SNP:

		if "*" in genos:
			return "pass"

		empty_set.add(genos[0])
		empty_set.add(genos[-1])


	if "." in empty_set: 
		empty_set.remove(".")

		print(empty_set)

	if len(empty_set) > 4:
		print( SNP ) 
		print( empty_set )
		print( "This pipeline can't really deal with that... so it's quitting. ")
		return "pass"
		# quit()

	marker_dict = {}
	for idx, elements in enumerate(empty_set):
		marker_dict[ elements ] = geno_list[idx]

	new_mkr_list = []
	for mkrs in SNP:

		new_mkr = ""
		for char in mkrs: 

			if (char == "/") or (char == "|"): 
				continue

			# if (char == ".") or (char == "*"):
			if (char == "*"): 
				new_mkr = "-"
				break
			
			else:
				new_mkr = new_mkr + marker_dict[char]

		new_mkr_list.append( new_mkr )

	return new_mkr_list  


def main(mod_count = 0, debug = False):

	"""This function converts "diploid" SNP data FROM 529L/P60002 parents and progeny 
		into 4way data for use in recombination_analysis.py """

	arg_list = []
	for arg in sys.argv:
		arg_list.append( arg )
	
	path_and_name = arg_list[1] 
	dir_path, file_name = parse_path( path_and_name )

	# CHECK THE PARENTS BEFORE YOU BEGIN!!
	parent_1 = "MAnderson017_1_R225-2_S31.GT"
	parent_2 = "MAnderson017_1454_R225-2_S101.GT"

	five29L_cross = {}

	header = 1

	with open( path_and_name, "r" ) as data_file:

		for idx, snps in enumerate(data_file): 

			snps = snps.strip().split("\t")
			
			if header: 
				header = 0
				
				parent_1_idx = snps.index(parent_1)
				parent_2_idx = snps.index(parent_2)
				
				print("Please confirm Parent 1 is index",parent_1_idx, "and Parent 2 is index", parent_2_idx)
				
				order_bool = (parent_1_idx < parent_2_idx)

				if order_bool:
					five29L_cross["header"] = ["Chromosome", "Chr_Position", "# Catalog ID"] + [snps[parent_1_idx]] + [snps[parent_2_idx]] + snps[9:parent_1_idx] + snps[parent_1_idx+1:parent_2_idx] + snps[parent_2_idx+1:] 
				else: 
					five29L_cross["header"] = ["Chromosome", "Chr_Position", "# Catalog ID"] + [snps[parent_1_idx]] + [snps[parent_2_idx]] + snps[9:parent_2_idx] + snps[parent_2_idx+1:parent_1_idx] + snps[parent_1_idx+1:]

				continue

			chromosome = snps[0]
			chr_pos    = snps[1]
			par_1      = snps[parent_1_idx]
			par_2      = snps[parent_2_idx]

			ref_allele = snps[3]
			alt_allele = snps[4]

			cat_ID = str(idx)

			# if "Ca19-mtDNA" in chromosome: # I need to include mitochondria
			# 	continue

			if order_bool: 
				prog_snps = [snps[parent_1_idx]] + [snps[parent_2_idx]] + snps[9:parent_1_idx] + snps[parent_1_idx+1:parent_2_idx] + snps[parent_2_idx+1:] 
			else: 
				prog_snps = [snps[parent_1_idx]] + [snps[parent_2_idx]] + snps[9:parent_2_idx] + snps[parent_2_idx+1:parent_1_idx] + snps[parent_1_idx+1:] 

			new_mkr = snp2mkr(prog_snps, cat_ID, ref_allele, alt_allele)

			mod_count += 1
			if (debug) and (mod_count % 1000 == 0):

				print(len(prog_snps), " Input: ", prog_snps)
				print(len(new_mkr), " Output: ", new_mkr)

			if new_mkr == "pass":
				continue

			new_529L_mkr = [chromosome, chr_pos, cat_ID] + new_mkr

			five29L_cross[idx] = new_529L_mkr 

	new_5_file = open( dir_path + "prog_4way.csv", "w" )

	csv_printer( five29L_cross, new_5_file )

main()
