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


def cf_graph_maker( allele_dict, file_path, cf_graph_prefix, file_names, only_one = 0 ):
		
	cf_graph_file = open( file_path + cf_graph_prefix + file_names.split("/")[-1].strip()[:-4]+".csv", "w" )
	for cat_nums, leel_info in allele_dict.items():
		for leels in leel_info: 
			print( ",".join([cat_nums]+leels+ [str(len(leel_info))] ), file = cf_graph_file )
			if only_one: 
				break

	cf_graph_file.close() 
	os.system( "Rscript cf_graph.R " + file_path + cf_graph_prefix + file_names.split("/")[-1].strip()[:-4]+".csv" )


def tags_file_to_dict( tags_file_name ):

	cat_tags_file = open( tags_file_name, "r" )
	cat_tags_list = file_reader( cat_tags_file, "\t" )

	cat_tag_dict = {}
	for lines in cat_tags_list: # Parse this file list into a dictionary
	    if "#" in lines[0]:
	        continue
	    new_key = lines[2]
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


def main( minimum_read_count = 100000, sort_output = 0, make_graphs = 0 ): 

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

	exclude_file = open( file_path + "excluded_strains.csv", "w" )
	list_file = open( file_path_n_name, "r" )

	for file_names in list_file: 
		
		leel_file = open( file_names.strip(), "r" )
		leel_dict = {}
		# Initate a different container for the ploidy investigations
		allele_dict = {}

		for lines in leel_file:
			if "#" in lines: 
				continue
			line_list = lines.strip().split("\t")

			if line_list[2] not in leel_dict: 
				leel_dict[ line_list[2] ] = int(line_list[-1].strip())
			else:
				leel_dict[ line_list[2] ] += int(line_list[-1].strip())

			if line_list[2] not in allele_dict: 
				allele_dict[ line_list[2] ] = []
				allele_dict[ line_list[2] ].append( line_list[3:] )

			else: 
				allele_dict[ line_list[2] ].append( line_list[3:] )

		# We don't do this, but calculating the number of different "aberrant" ( 3+ alleles in a marker ) could be beneficial 
		if make_graphs:
			cf_graph_maker( allele_dict, file_path, cf_graph_prefix, file_names, only_one = 0 )
		
		# Combine the data an make a new directory to store all the new files in 
		mismatches = 0
		prog_dict = {}
		prog_list = []
		for cat_nums, tags in tags_dict.items(): 
			
			if (cat_nums not in prog_dict) and (cat_nums in leel_dict): 
				prog_dict[ cat_nums ] = tags_dict[ cat_nums ] + [str(leel_dict[ cat_nums ])]
				prog_list.append( [cat_nums] + [tags_dict[ cat_nums ][0]] + [int( tags_dict[ cat_nums ][1] )] + [str(leel_dict[ cat_nums ])] )
			
			# Add lines for empty markers 
			elif (cat_nums not in prog_dict) and (cat_nums not in leel_dict): 
				prog_dict[ cat_nums ] = tags_dict[ cat_nums ] + [ "0" ]
				prog_list.append( [cat_nums] + [tags_dict[ cat_nums ][0]] + [int( tags_dict[ cat_nums ][1] )] + ["0"] ) # Turn to integers for sorting

			# Make a list, sort them by chromosomal location. 
			if sort_output: 
				prog_list = sorted( prog_list, key=operator.itemgetter(1,2) )

		new_file = open( file_path + new_file_prefix + file_names.split("/")[-1].strip()[:-4]+".csv", "w" )
		
		read_sum = 0

		for lines in prog_list:
			str_line = []
			read_sum += int(lines[-1])
			for elements in lines: 
				str_line.append( str(elements) ) # Turn the ints back into strings so they can be joined and printed easily. 
			printable_line = ",".join(str_line).strip()
			print( printable_line, file = new_file )

		if read_sum < minimum_read_count: 
			print( file_names.strip(), file = exclude_file )

		new_file.close()
		leel_file.close()

	# Simplest solution here is to run the ploidy files into R to get the sums (the data is in a single column). If the file's too big, 
	os.system( "rm " + file_path_n_name ) 
	exclude_file.close()

main( sort_output = 0 )

# Command: 
# python3 allele_analysis.py SCxP60002/batch_2.catalog.tags.tsv SCxP60002/allele_files_for_ploidy_measurements/