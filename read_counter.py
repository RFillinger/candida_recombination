
import subprocess
import sys
import os

def main(): 

	arg_list = []
	for arg in sys.argv:
		arg_list.append( arg )

	dir_path = arg_list[1]

	if dir_path[-1].strip() != "/": 
		dir_path = dir_path.strip() + "/"

	# Make a file with all the names of the .bam files in the directory of execution
	os.system( "ls " + dir_path + "*.bam* > " + dir_path + "temporary.txt" )

	# Open the file with everything in it to read all the lines in it
	file_list = open( dir_path + "temporary.txt", "r" )
	read_count_file = open( dir_path + "read_counts.txt", "w" )

	# This is for parsing readcounts from strings later
	numerals = "0123456789"

	for files in file_list: 

		read_count_obj = subprocess.run(["samtools", "view", "-c", "-F", "260", files], stdout=subprocess.PIPE)
		read_count_str = read_count_obj.stdout

		simple_read_cnt_str = ""
		for char in str(read_count_str): 
			if char in numerals: 
				simple_read_cnt_str += char

		print( files.split(".")[0] + "," + simple_read_cnt_str, file = read_count_file )
		print( files.split(".")[0] + "," + simple_read_cnt_str )

	os.system( "rm " + dir_path + "temporary.txt" )

	file_list.close()
	read_count_file.close()

main()