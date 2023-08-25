import time
import sys

def file_reader( input_file, sep = "," ):

	'''Reads in .csv files into a list'''

	line_list = []
	for lines in input_file: 
		line_list.append( lines.strip().split(sep) )

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

	# Convert the data to 4way here
	geno_list = ["a", "b", "c", "d"]

	empty_set = set()
	for genos in SNP[0:2]: # Markers from the parents

		if "/" in genos: 
			genotype_list = genos.split("/")
		elif "|" in genos: 
			genotype_list = genos.split("|")

		empty_set.add(genotype_list[0])
		empty_set.add(genotype_list[1])


	if len(empty_set) > 4:
		# print(SNP) 
		# print( empty_set )
		# print( "There are 5 alleles at the marker:", chromosome, chr_pos)
		# print( "This pipeline can't really deal with that... so it's quitting. ")
		return ""

	marker_dict = {}
	for idx, elements in enumerate(empty_set):
		marker_dict[ elements ] = geno_list[idx]	

	new_mkr_list = []
	for mkrs in SNP:

		new_mkr = ""

		for genotypes in mkrs:

			if (genotypes == "/") or (genotypes == "|"): 
				continue

			try: 
				new_mkr += marker_dict[ genotypes ]

			except KeyError: 

				new_mkr = "-"
				break					

		new_mkr_list.append( new_mkr )

	return new_mkr_list  


def main():

	"""This function converts "diploid" SNP data FROM 529L/P60002 parents and progeny 
		into 4way data for use in recombination_analysis.py """

	arg_list = []
	for arg in sys.argv:
		arg_list.append( arg )
	
	path_and_name = arg_list[1] 
	dir_path, file_name = parse_path( path_and_name )
		
	data_file = open( file_name, "r" )
	data_list = file_reader( data_file, sep = "\t" )

	cross = []
	par1_col_idx = 28
	par2_col_idx = 23

	header = 1
	for index, snps in enumerate(data_list): 

		try: 
			chromosome        = snps[0]
			chr_pos           = snps[1]
			cat_ID = str(index)

		except IndexError: 
			print(data_list[index])

		if  ("Ca19-mtDNA" in chromosome): 
			continue

		if header: 
			header = 0
			header_line = ",".join(["Chromosome", "Chr_Position", "# Catalog ID", snps[par1_col_idx], snps[par2_col_idx]] + \
				snps[9:min(par1_col_idx, par2_col_idx)] + \
				snps[min(par1_col_idx, par2_col_idx)+1:max(par1_col_idx, par2_col_idx)] + \
				snps[max(par1_col_idx, par2_col_idx)+1:])
			cross.append( ["Chromosome", "Chr_Position", "# Catalog ID", snps[par1_col_idx], snps[par2_col_idx]] + \
				snps[9:min(par1_col_idx, par2_col_idx)] + \
				snps[min(par1_col_idx, par2_col_idx)+1:max(par1_col_idx, par2_col_idx)] + \
				snps[max(par1_col_idx, par2_col_idx)+1:] )

			print( "Parent 1: ", snps[par1_col_idx], "\nParent 2: ", snps[par2_col_idx])
			continue

		SCx529L = [ snps[par1_col_idx], snps[par2_col_idx] ] + \
				snps[9:min(par1_col_idx, par2_col_idx)] + \
				snps[min(par1_col_idx, par2_col_idx)+1:max(par1_col_idx, par2_col_idx)] + \
				snps[max(par1_col_idx, par2_col_idx)+1:] # Puts the parents in the beginning and keeps remaining progeny

		new_mkr = []

		marker = snp2mkr(SCx529L, cat_ID)
		
		if marker == "": # returned nothing due to error: 
			continue
		
		new_mkr.append(marker)

		new_529L_mkr = [chromosome, chr_pos, cat_ID] + new_mkr[0]

		cross.append( new_529L_mkr )

	new_5_file = open( "all_prog_WGS.csv", "w" )
	csv_printer( transpose(cross), new_5_file )


main()
