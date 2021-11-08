import os

# renames files like P60003 and SC5315 back to SC5314. Sorta meant for a one time use.

sc_file = open("sc_change.txt", "r")
p6_file = open( "p6_change.txt", "r")

for file_names in sc_file: 
	line_list = file_names.split("_")
	line_list[0] = "SC5314"
	name_string = "_".join(line_list)
	print( "mv " + file_names.strip() + " " + name_string )
	os.system("mv " + file_names.strip() + " " + name_string)

for file_names in p6_file: 
	line_list = file_names.split("_")
	line_list[0] = "P60002"
	name_string = "_".join(line_list)
	print( "mv " + file_names.strip() + " " + name_string)
	os.system("mv " + file_names.strip() + " " + name_string)