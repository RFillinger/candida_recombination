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

	    mkr_loc = lines[3] + ":" + lines[4]
	    cat_tag_dict[mkr_loc] = lines[2] # Stores marker catalog numbers for their respective locations; just to easily input them into the new file

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

	if file_path[-1] == "/":
		temp_file_name = "temporary_1.txt"
		file_path_n_name = file_path + temp_file_name
		
		temp_file_name_2 = "temporary_2.txt"
		file_path_n_name_2 = file_path + temp_file_name_2
		
	else:
		temp_file_name = "temporary_1.txt"
		file_path_n_name = file_path + "/" + temp_file_name
		
		temp_file_name_2 = "temporary_2.txt"
		file_path_n_name_2 = file_path + "/" + temp_file_name_2
		
	file_list_preparer( file_path, ending = "*.alleles.tsv", out_file_name = temp_file_name  )
	file_list_preparer( file_path, ending = "*.tags.tsv", out_file_name = temp_file_name_2  )

	if exclude_chr: 
		exclude_file = open( "progeny_ploidy_blacklist.csv", "r" ) 
		exclude_dict = {}
		for lines in exclude_file: 
			line_list = lines.split(",")
			samples, chromosome_str = line_list[0], line_list[1].strip()
			chromosome_list = chromosome_str.split(",")
			exclude_dict[ samples ] = chromosome_list

	# Get the data for moving location data into final markers to remove them more effeciently
	tags_dict = tags_file_to_dict( tags_file_name )

	# Pull in information from "*.alleles.tsv"
	list_file = open( file_path_n_name, "r" )

	progeny_dict = {}

	for file_names in list_file: 
		sample_name = file_names.split( "/" )[-1][:-13].strip()
		progeny_dict[ sample_name ] = {}
		leel_file = open( file_names.strip(), "r" )

		# Initate a different container for the ploidy investigations
		for lines in leel_file:
			if "#" in lines: 
				continue
			line_list = lines.strip().split("\t")
			line_list[4] = float(line_list[4])
			catalog_num = int(line_list[2].strip())
			if catalog_num not in progeny_dict[ sample_name ]: 
				progeny_dict[ sample_name ][ catalog_num ] = {"allele":[], "position":[]}
				progeny_dict[ sample_name ][ catalog_num ]["allele"].append( line_list[3:]) 
			else: 
				progeny_dict[ sample_name ][ catalog_num ]["allele"].append( line_list[3:])
					
		leel_file.close()
	list_file.close()

	list_file_2 = open( file_path_n_name_2, "r" )

	# Enter information for the locations for preparations for the second file
	for file_names in list_file_2: 
		sample_name = file_names.split( "/" )[-1][:-10].strip()
		if "batch" in sample_name: 
			continue
		leel_file = open( file_names.strip(), "r" )

		matches = 0
		total = 0

		# Initate a different container for the ploidy investigations
		for lines in leel_file:

			if "#" in lines: 
				continue
			
			line_list = lines.strip().split("\t")
			if "Ca21chr" not in line_list[3]: 
				continue

			catalog_num = int(line_list[2].strip())
			chromosome = line_list[3]
			pos_str = line_list[4]

			position = [chromosome, pos_str]

			if catalog_num in progeny_dict[ sample_name ]:
				matches += 1
				total += 1
				progeny_dict[ sample_name ][ catalog_num ]["position"] = position
			else: 
				total += 1

		# print( "\nSample: ", sample_name )
		# print( "Total Markers: ", total )
		# print( "Het Marker Matches: ", matches )
				
		leel_file.close()

	list_file_2.close()

	header = "sample,catalog_ID,chromosome,chrom_loc,total_leel_count,total_depth,num_1_leel_prop,num_2_leel_prop"
	if exclude_chr: 
		new_file = open( file_path + "excluded_markers.csv", "w" )

	else: 
		new_file = open( file_path + "all_inclusive_markers.csv", "w" )
	print( header, file = new_file )

	for sample_name in progeny_dict:

		if exclude_chr:
			exclusionary_chromosomes = []
			if (sample_name in exclude_dict):
				if exclude_dict[ sample_name ][0] == "": 
					continue
				else: 
					for chromosomes in exclude_dict[ sample_name ]:
						new_chromosome = "Ca21chr" + chromosomes.replace("\"", "") + "_C_albicans_SC5314" # This would have to be fixed in the 
						exclusionary_chromosomes.append( new_chromosome )

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
				chrom = info_dict["position" ][0]
				chrom_pos = info_dict["position"][1]

				tags_dict_key = chrom + ":" + chrom_pos
				if tags_dict_key in tags_dict: 
					univ_tags_mkr = tags_dict[ tags_dict_key ]

				else: 
					univ_tags_mkr = "missing"

				line_list = [ sample_name, univ_tags_mkr, chrom, chrom_pos, str(allele_count), str(total_read_depth) ] + allele_read_prop_list
				new_line = ",".join( line_list )

				if exclude_chr: 
					if chrom in exclusionary_chromosomes: 
						continue
				
				# if (total_read_depth > 0) and ( univ_tags_mkr != "missing"): 
				if (total_read_depth > 0):
					print(new_line, file = new_file )

			except IndexError: # This screens for heterozygous markers; markers with only one allele
				pass

	new_file.close()

	# Simplest solution here is to run the ploidy files into R to get the sums (the data is in a single column). If the file's too big, 
	os.system( "rm " + file_path_n_name + " " + file_path_n_name_2)

start = 0
while start < 2: 
	main( exclude_chr = start )
	start += 1

# Sample commands: 
# python3 allele_demystifier.py SCx529L/batch_1.catalog.tags.tsv SCx529L/allele_files_for_ploidy_measurements/ 
# python3 allele_demystifier.py SCxP60002/batch_2.catalog.tags.tsv SCxP60002/allele_files_for_ploidy_measurements/