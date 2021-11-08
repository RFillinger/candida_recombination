import sys
import os
import statistics
import random

# Somethings that need to be done: 
# 1)  Sort the output transposed genotypes pile by chromosome and then by chromosome postion so you don't have to do it in excel
# 2)  Take the .catalog.tags.tsv file in and transpose it yourself using numpy so that there are fewer complications
# 3)  Replace the names of the chromosomes with numbers of the chromosomes
# 4)  Integrate the phenotypes into the output file based off of input
# 5)  Replace "--" with "-" and uppercase "AB" alleles with lowercase "ab" alleles
# 6)  Combine all the genotypes data analysis tools into a single command
# 7)  Prepare an input file
# 8)  Clean out functions that you've deemed outdated and/or useless, at least comment them out or label them as such
# 9)  Walk your code and document it, extend your markdown file to accomodate any changes in the command for this program
# 10) Statistics on the number of markers per progeny and perhaps an output file to analyze the reads in R (or output the graph yourself)
# 11) Find the average number of samples that all markers occupy (that's the number of progeny you can use for the marker screen)
# 12) FIX THE F2 CONVERSION!!!!!! IT DOESN'T WORK!!!!!!

arg_list = []

for arg in sys.argv:
    arg_list.append( arg )

input_file_name = arg_list[1] #File of transposed .genotypes file from STACKS
tags_file_name  = arg_list[2] #Catalog.tags file that has been transposed from the same stacks output
parent_1_str    = arg_list[3] #String that is found exclusively in parent 1
parent_2_str    = arg_list[4] #String found exclusively in parent 2
output_name     = arg_list[5] #Name of the output
#view_for_recombination = arg_list[6] #Whether or not you want the data to be interpretable for recombination (only alleles that are different between parents)


def tags( markers, tags_array ):

    index_list = []
    index_list.append(0)

    magnus_tags_array = []

    count = 0
    for geno_markers in markers: 
        try:
            index_list.append( tags_array[0].index( geno_markers ) )
        except ValueError:
            #Line ID will trip this
            continue
    
    tags_mrkr_list    = []
    tags_chr_list     = []
    tags_chr_pos_list = []
    tags_strand_list  = []
    tags_seq_list     = []

    for indices in index_list:

        tags_mrkr_list   .append( tags_array[0][ indices ] )
        tags_chr_list    .append( tags_array[1][ indices ] )
        tags_chr_pos_list.append( tags_array[2][ indices ] )
        tags_strand_list .append( tags_array[3][ indices ] )
        tags_seq_list    .append( tags_array[4][ indices ] )

    magnus_tags_array.append(  tags_mrkr_list   )
    magnus_tags_array.append(   tags_chr_list   )
    magnus_tags_array.append( tags_chr_pos_list )
    magnus_tags_array.append( tags_strand_list  )
    magnus_tags_array.append(   tags_seq_list   )

    return magnus_tags_array


def multi_find( main_string, str_to_find ):

    main_string_length = len(main_string)

    split_size = len(str_to_find)

    i = 0

    index_list = []

    while (i + split_size) <= main_string_length:

        main_string_split = main_string[i:i+split_size]

        if str_to_find in main_string_split:

            index_list.append(i)

        i += 1

    return index_list


def binary_progeny_alleles( progeny, diff_alleles ):

    '''Extracts indexes from "progeny" list given by "diff_alleles" 
    (different alleles which are alleles in the progeny that can be traced to specific parents)'''
    
    new_array = []
    progeny_match_list = []
    
    progeny_match_list.append( progeny[0] )

    for indices in diff_alleles:

        progeny_match_list.append( progeny[ indices ] )

    return progeny_match_list


def parent_matcher( file_array ):

    parent_1 = file_array[0].split( "," )
    parent_2 = file_array[1].split( "," )

    out_list = []

    for progs in file_array[2:]:

        progs = progs.split(",")
        new_prog = ""
        index = 0
        while index < len( parent_1 ):

            if index == 0:
                new_prog = progs[index] + "," #add the sample name
                index += 1
                continue

            if "-" in progs[index].lower():
                new_prog += "-,"
                index +=1
                continue

            if progs[index].lower() == parent_1[index].lower():
                new_prog += "1,"
                index += 1

            elif progs[index].lower() == parent_2[index].lower():
                new_prog += "2,"
                index += 1

            else:
                new_prog += "n,"
                index += 1

        else:
            new_prog = new_prog[:-1] #removes last comma for csv

        out_list.append( new_prog )

    return out_list


def find_het_indices( list_1, list_2, start_index = 0, absolute_len = 2, ignore_val = "XX" ):

    '''Finds matching or partially matching indexes shared by two equally sized lists

        absolute_len:  is the ploidy for the organism, default is 2 (C. albicans is diploid)
        ignore_val:    any alleles that you want to be ignored by the function'''

    matching_indices = []
    len_match_bool = len( list_1 ) == len( list_2 )

    if len_match_bool: 

        list_length = len( list_1 ) - 1

        index = start_index

        while index <= list_length:

            value_1 = list_1[ index ].lower()
            value_2 = list_2[ index ].lower()

            if ( len( value_1 ) != absolute_len ) or \
               ( len( value_2 ) != absolute_len ):
                index += 1 
                continue

            if ( value_1 == ignore_val ) or \
               ( value_2 == ignore_val ):
                index += 1 
                continue

            empty_set = set()

            empty_set.add( value_1[0] )
            empty_set.add( value_1[1] )
            empty_set.add( value_2[0] )
            empty_set.add( value_2[1] )

            if len( empty_set ) >= 3: #If there are 3 unique values, keep it
                matching_indices.append( index )
                index += 1 
                continue

            else: 

                if ( value_1[0] == value_1[1] ) and \
                   ( value_2[0] == value_2[1] ) and \
                   ( value_1[0] != value_2[0] ):
                    matching_indices.append( index )
                    index += 1 
                    continue

                else:
                  index += 1 
                  continue                        

            index += 1

        return matching_indices

    else: 

        print( "ERROR in find_het_indices: \n" \
               "LIST SIZES ARE NOT EQUAL")
        quit()

#######################################################################################
### Input file preparation for easier visualization

input_file = open( input_file_name, "r" )

file_array = []
parent_1   = []
parent_2   = []

# Read the file into 3 arrays (two seperate parent arrays and an array of all progeny) that can be used
# The two parent arrays make conversions of them easier and also simplifies the readability a bit 
for lines in input_file: 

    prepped_line = lines.strip().split(",")
    
    if multi_find( prepped_line[0], parent_1_str ) != []:
        parent_1.append( prepped_line )

    elif multi_find( prepped_line[0], parent_2_str ) != []:
        parent_2.append( prepped_line )

    else: #Removes parents from the list so only progeny are used
        file_array.append( prepped_line )

input_file.close()

#######################################################################################

#This next chunk of code finds the most common allele for each parent at each position
marker_start_idx = 0 
sample_start_idx = 6

num_markers = len( file_array[0] ) - 1 #Total number of markers in the analysis
num_samples = len(  file_array   ) - 6 #Total number of samples

p1_error_indexes = [] # I don't use these beyond data collection, but they may be important in future stacks studies
p2_error_indexes = []

parent_1_mode = []
parent_2_mode = []

parent_1_cat = []
parent_2_cat = []

parent_1_cat.append( "Parent_1_" + parent_1_str )
parent_2_cat.append( "Parent_2_" + parent_2_str )

while marker_start_idx <= num_markers: 

    temp_list_1 = []
    temp_list_2 = []
    
    for parents in parent_1:

        temp_list_1.append( parents[ marker_start_idx ] )
        
    for other_parents in parent_2:

        temp_list_2.append( other_parents[ marker_start_idx ] )


    if marker_start_idx > 0:

        ran_bool = random.randint( 0, 1 ) # For picking the mode where one des not exist (not yet implemented; not simple)

        try:

            p1_marker_mode = statistics.mode( temp_list_1 )
            parent_1_cat.append( p1_marker_mode )

        except statistics.StatisticsError: 

            p1_error_indexes.append( marker_start_idx )
            parent_1_cat.append( "-" )
        
        try:

            p2_marker_mode = statistics.mode( temp_list_2 )
            parent_2_cat.append( p2_marker_mode )

        except statistics.StatisticsError: 
            p2_error_indexes.append( marker_start_idx )
            parent_2_cat.append( "-" ) #There are very few of these 
    
    marker_start_idx += 1

#######################################################################################

tags_file = open( tags_file_name, "r" ) #opens catalog file containing marker information

empty_array = []
for lines in tags_file:

    tags_line = lines.split(",")
    empty_array.append( tags_line )

tags_array = empty_array[2:6]
tags_array.append( empty_array[9] )

# for lines in tags_array: 
#     print( lines[0:10] )

tags_array[1].insert( 0, "Chromosome"   )
tags_array[2].insert( 0, "Chr_Position" )
tags_array[0].insert( 0, "Marker_Num"   )
tags_array[3].insert( 0, "Strandedness" )
tags_array[4].insert( 0, "Sequence"     )

new_file = open( output_name, "w" )

tags_array = tags( file_array[0], tags_array )

tags_chr     = tags_array[1]
tags_chr_pos = tags_array[2]
tags_markers = tags_array[0] 
tags_strand  = tags_array[3] # This shows strandedness in the output. Not really necessary, but here you go.
tags_seq     = tags_array[4]

lines_array = []

lines_array.append( ",".join( tags_chr      ))
lines_array.append( ",".join( tags_chr_pos  ))
lines_array.append( ",".join( file_array[0] ))
# lines_array.append( ",".join( tags_strand   ))
lines_array.append( ",".join( parent_1_cat  ))
lines_array.append( ",".join( parent_2_cat  ))

for progeny in file_array[6:]:

    lines_array.append( ",".join( progeny ))

for elements in lines_array:
    print( elements, file = new_file )

for addons in parent_matcher( lines_array[4:] ):
    print( addons, file = new_file )

new_file.close()
tags_file.close()