
from subprocess import Popen, PIPE, CalledProcessError
import sys
import os

def main(): 

	arg_list = []
	for arg in sys.argv:
		arg_list.append( arg )

	dir_path = arg_list[1]
	input_bed = arg_list[2]

	if dir_path[-1].strip() != "/": 
		dir_path = dir_path.strip() + "/"

	# Make a file with all the names of the .bam files in the directory of execution
	os.system( "ls " + dir_path + "*.bam* > temporary.txt" )

	# Open the file with everything in it to read all the lines in it
	temp_file = open( "temporary.txt", "r" )
	read_count_file = open( dir_path + "read_counts.txt", "w" )

	# This is for measuring progress
	read_list_number_file = open( "temporary.txt", "r" )
	number = len( read_list_number_file.readlines() )
	
	print()

	# Define the dictionary that will contain all our data
	read_count_dictionary = {}
	read_count_dictionary[ "header" ] = []
	# For testing: 
	test = 0

	if test: 
		temp_file = open( "perma_temp.txt", "r" )

	header_list = []
	file_count = 1

	for files in temp_file:

		# Notify user of progress: 
		print( "Counting read depths in " + files.strip() + " (" + str(file_count) + "/" + str(number) + ")" )
		file_count += 1

		cmd = [ "bedtools","coverage", "-counts", "-a", input_bed, "-b", files.strip() ]
		# Printing the command to check if testing:  
		if test:
			print( " ".join( cmd ))

		sample_name = files.split("/")[-1].split(".")[0]
		read_count_dictionary[ "header" ].append( sample_name )

		with Popen( cmd, stdout=PIPE, bufsize = 1, universal_newlines = True ) as p:
			for line in p.stdout:
				line_list = str(line).strip().split("\t")
				marker_key = line_list[0] + ":" + line_list[1] + ":" + line_list[2] 
				feature_depth = line_list[6]

				if marker_key not in read_count_dictionary:
					read_count_dictionary[ marker_key ] = [ feature_depth.strip() ]

				# Check to make sure the size iw what you want.  
				elif marker_key in read_count_dictionary: 
					if len( read_count_dictionary[ "header" ]) == len(read_count_dictionary[marker_key])+1: 
						read_count_dictionary[ marker_key].append( feature_depth.strip() )


	for markers, read_counts in read_count_dictionary.items():

		if markers == "header": 
			print( str( "Chromosome\tMkrStart\tMkrEnd\t") + "\t".join(read_counts), file = read_count_file ) 
			continue
		print( str( "\t".join(markers.split(":")) + "\t" + "\t".join(read_counts) ), file = read_count_file )

	os.system( "rm temporary.txt" )

	temp_file.close()
	read_count_file.close()


main()