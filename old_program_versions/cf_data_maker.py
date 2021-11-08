# The goal of this program is to take individual allele files from each strain and combine them with the catalog tags file. Then, I can give the output to Anna
import sys
import os
import operator

def file_reader( input_file, delimiter = "," ):
    '''Reads in a file to a list'''
    line_list = []
    for lines in input_file: 
        line_list.append( lines.strip().split(delimiter) )
    return line_list


def tags_file_to_dict( tags_file_name ):

	cat_tags_file = open( tags_file_name, "r" )
	cat_tags_list = file_reader( cat_tags_file, "\t" )

	cat_tag_dict = {}
	for lines in cat_tags_list: # Parse this file list into a dictionary
	    if "#" in lines[0]:
	        continue
	    new_key = lines[2].strip()
	    cat_tag_dict[new_key] = lines[3:5]

	return cat_tag_dict


# First function I need to make will report all the file allele file names and put them into a file that this program can read
def file_list_preparer( file_path, ending = "*.csv", out_file_name = "temporary.txt" ):

	""" Puts all the names of files containing ending into a text file """
	if file_path[-1] == "/":
		file_path_n_name = file_path + out_file_name
		os.system( "ls " + file_path + ending + " > " + file_path_n_name )
	else:
		file_path_n_name = file_path + "/" + out_file_name
		os.system( "ls " + file_path + "/" + ending + " > " + file_path_n_name )


def main( exclude_chr = 1 ): 

	arg_list = []
	for arg in sys.argv:
		arg_list.append( arg )

	tags_file_name = arg_list[1]
	file_path = arg_list[2]

	new_file_prefix = "ploidy_"
	cf_graph_prefix = "cf_"

	tags_dict = tags_file_to_dict( tags_file_name )

	temp_file_name = "temporary.txt"
	file_list_preparer( file_path, ending = "*.alleles.tsv", out_file_name = temp_file_name  )
	if file_path[-1] == "/":
		file_path_n_name = file_path + temp_file_name
	else:
		file_path_n_name = file_path + "/" + temp_file_name

	list_file = open( file_path_n_name, "r" )

	if exclude_chr: 
		exclude_file = open( "progeny_ploidy_blacklist.csv", "r" )
	else: 
		exclude_file = open( "blank_blacklist.csv", "r" )

	if exclude_chr: 
		exclude_dict = {}
		for lines in exclude_file: 
			line_list = lines.split(",")
			samples, chromosome_str = line_list[0], line_list[1].strip()
			chromosome_list = chromosome_str.split(",")
			exclude_dict[ samples ] = chromosome_list

	exclude_file.close()

	progeny_dict = {}

	for file_names in list_file: 
		sample_name = file_names.split( "/" )[-1][:-13]

		if exclude_chr:
			exclusionary_chromosomes = []
			if (sample_name in exclude_dict):
				if exclude_dict[ sample_name ][0] == "": 
					continue
				else: 
					for chromosomes in exclude_dict[ sample_name ]:
						new_chromosome = "Ca21chr" + chromosomes.replace("\"", "") + "_C_albicans_SC5314"
						exclusionary_chromosomes.append( new_chromosome )

		progeny_dict[ sample_name ] = {}

		leel_file = open( file_names.strip(), "r" )
		leel_dict = {}
		# Initate a different container for the ploidy investigations
		allele_dict = {}

		for lines in leel_file:

			if "#" in lines: 
				continue

			line_list = lines.strip().split("\t")
			line_list[4] = float(line_list[4])

			catalog_num = line_list[2] 
			if catalog_num not in progeny_dict[ sample_name ]: 
				progeny_dict[ sample_name ][ catalog_num ] = {"allele":[], "position":[]}
				progeny_dict[ sample_name ][ catalog_num ]["allele"].append( line_list[3:]) 
			else: 
				progeny_dict[ sample_name ][ catalog_num ]["allele"].append( line_list[3:]) 
		
		for cat_nums, tags in tags_dict.items(): 
			chromosome = tags[0]
			for sample_name, mkr_catalog_dict in progeny_dict.items(): 
				
				if exclude_chr: 
					if chromosome in exclusionary_chromosomes: 
						progeny_dict[ sample_name][cat_nums] = {"allele":[], "position":[]}
						progeny_dict[ sample_name][cat_nums]["allele"] = []
						progeny_dict[ sample_name][cat_nums]["position"] = tags

					else: 
						try: 
							progeny_dict[sample_name][ cat_nums ][ "position" ] = tags
						except KeyError: 
							progeny_dict[ sample_name][cat_nums] = {"allele":[], "position":[]}
							progeny_dict[ sample_name][cat_nums]["allele"] = []
							progeny_dict[ sample_name][cat_nums]["position"] = tags
				else: 
					try: 
						progeny_dict[sample_name][ cat_nums ][ "position" ] = tags
					except KeyError: 
						progeny_dict[ sample_name][cat_nums] = {"allele":[], "position":[]}
						progeny_dict[ sample_name][cat_nums]["allele"] = []
						progeny_dict[ sample_name][cat_nums]["position"] = tags

		leel_file.close()

	if exclude_chr: 
		new_file = open( file_path + "excluded_markers.csv", "w" )
	else: 
		new_file = open( file_path + "all_inclusive_markers.csv", "w" )

	for sample_name in progeny_dict:
		for cat_nums, info_dict in progeny_dict[ sample_name ].items(): 
			allele_count = len( progeny_dict[ sample_name ][ cat_nums ][ "allele" ] )
			srt_alleles = sorted( progeny_dict[ sample_name ][ cat_nums ][ "allele" ], key=operator.itemgetter(1), reverse = 1 )
			# srt_alleles = progeny_dict[ sample_name ][ cat_nums ][ "allele" ] # No sorting
			total_read_depth = 0
			two_stop = 0
			allele_read_prop_list = []
			for alleles in srt_alleles: 
				total_read_depth += int( alleles[2] )
				# I only want the top 2
				if two_stop < 2: 
					allele_read_prop_list.append( str(alleles[1]) )
					two_stop += 1

			try: 
				chrom = info_dict[ "position" ][0]
				chrom_pos = info_dict[ "position"][1]
				line_list = [ sample_name, cat_nums, chrom, chrom_pos, str(allele_count), str(total_read_depth) ] + allele_read_prop_list
				new_line = ",".join( line_list )
			except IndexError: 
				# line_list = [ sample_name, cat_nums, "NA", "NA", str(total_read_depth) ] + allele_read_prop_list 
				# new_line = ",".join( line_list )
				pass
			
			if total_read_depth > 0: 
				print(new_line, file = new_file )


	new_file.close()

	# Simplest solution here is to run the ploidy files into R to get the sums (the data is in a single column). If the file's too big, 
	os.system( "rm " + file_path_n_name ) 

start = 0
while start < 2: 
	main( exclude_chr = start )
	start += 1


# Command: 
# python3 cf_data_maker.py SCxP60002/batch_2.catalog.tags.tsv SCxP60002/allele_files_for_ploidy_measurements/