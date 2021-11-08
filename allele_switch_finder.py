import sys
import os
import statistics
import random

def file_reader( input_file, delimiter = "," ):
    '''Reads in a file to a list'''
    line_list = []
    for lines in input_file: 
        line_list.append( lines.strip().split(delimiter) )
    return line_list


def csv_printer( csv_list, output_file ):

    for items in csv_list: 
        line = ",".join(items)
        print( line, file = output_file)


def transpose(l1):

    new_list = []
    # iterate over list l1 to the length of an item
    for i in range(len(l1[0])):
        # print(i)
        row =[]
        for item in l1:
            # appending to new list with values and index positions
            # i contains index position and item contains values
            row.append(item[i])
        new_list.append(row)
    return new_list


def parent_combiner( par_prog_dict, par_1_str, par_2_str, header_key = "# Catalog ID" ): 

    p1_indices = []
    p2_indices = []
    prog_indices = []
    temp_dict = {}
    p1_dict = {}
    p2_dict = {}

    for i, names in enumerate(par_prog_dict[ header_key ]): 
        if par_1_str in names:
            p1_indices.append( i )
        elif par_2_str in names: 
            p2_indices.append( i )
        else: 
            prog_indices.append( i ) # The goal here is to separate the parents, select the most frequent allele from the parents and then filter out the progeny

    # print( "Parent 1: ", p1_indices )
    # print( "Parent 2: ", p2_indices )
    # print( "Progeny: ", prog_indices )

    for prog, mkrs in par_prog_dict.items(): 
        new_list = []
        p1_list  = []
        p2_list  = []
        for index, mkr in enumerate(mkrs):
            if index in prog_indices: 
                new_list.append( mkr ) # This will include the header information and 
            elif index in p1_indices: 
                p1_list.append( mkr )
            else: 
                p2_list.append( mkr )
        temp_dict[ prog ] = new_list
        p1_dict[ prog ]  = p1_list 
        p2_dict[ prog ]  = p2_list

    # Next thing to do is combine the parents together into a consensus
    p1_no_modes = {} # These are for storing the markers with multiple options, so they can be proofed if needed
    p2_no_modes = {}

    mode_p1_dict = {}
    for keys, mkrs in p1_dict.items(): 
        mkr_list = [i for i in mkrs if i != "-"]
        if len( mkr_list ) == 0: # If the marker is empty, fill it with emptiness and keep on moving! 
            mode_p1_dict[ keys ] = "-"
            continue
        try:
            mode_p1_dict[ keys ] = statistics.mode(mkr_list) # The new consensus marker is most common of all the parent markers
        except statistics.StatisticsError:
            p1_no_modes[ keys ] = mkr_list 
            # mode_p1_dict[ keys ] = p1_dict[ keys ][ random.randint(0, (len(mkr_list)-1)) ] # If there is no mode, it's chosen randomly from the possible ones (probably a bad idea...)
            mode_p1_dict[ keys ] = "-" # Just remove it

    mode_p2_dict = {}
    for keys, mkrs in p2_dict.items():
        mkr_list = [i for i in mkrs if i != "-"]
        if len( mkr_list ) == 0: 
            mode_p2_dict[ keys ] = "-"
            continue
        try:
            mode_p2_dict[ keys ] = statistics.mode(mkr_list)
        except statistics.StatisticsError: 
            p2_no_modes[ keys ] = mkr_list
            # mode_p2_dict[ keys ] = p2_dict[ keys ][ random.randint(0, (len(mkr_list)-1)) ] 
            mode_p2_dict[ keys ] = "-" 

    # Now, reattach the parents to the original output file
    new_dict = {}
    for mkr_num, mkr_stack in temp_dict.items(): 

        if mkr_num == "# Catalog ID": 
            new_dict[ mkr_num ] = temp_dict[ mkr_num ][0:3] + [par_1_str] + [par_2_str] + temp_dict[ mkr_num ][3:]
            continue

        new_dict[ mkr_num ] = temp_dict[ mkr_num ][0:3] + [mode_p1_dict[ mkr_num ]] + [mode_p2_dict[ mkr_num ]] + temp_dict[ mkr_num ][3:]

    print( "There are ", len( p1_no_modes ), " inconclusive markers for Parent 1.") # Report the number of "inconclusive" markers from the parents
    print( "There are ", len( p2_no_modes ), " inconclusive markers for Parent 2.")

    return new_dict


def main( shoot_trouble = 0 ):

    if not shoot_trouble: 
        arg_list = []
        for arg in sys.argv:
            arg_list.append( arg )
        geno_file_name = arg_list[1] # .genotypes file from STACKS
        tags_file_name = arg_list[2] # Catalog.tags file that has been transposed from the same stacks output
        parent_1_str   = arg_list[3] # String that is found exclusively in parent 1
        parent_2_str   = arg_list[4] # String found exclusively in parent 2
        output_name    = arg_list[5] # Name of the output

    elif shoot_trouble == 1:
        geno_file_name = "SCxP60002/batch_2.genotypes_1.tsv"
        tags_file_name = "SCxP60002/batch_2.catalog.tags.tsv"
        parent_1_str   = "SC5314"
        parent_2_str   = "P60002"
        output_name    = "SCxP60002/P60002xSC_four_way.csv"

    elif shoot_trouble == 2:
        geno_file_name = "SCx529L/batch_1.genotypes_1.tsv"
        tags_file_name = "SCx529L/batch_1.catalog.tags.tsv"
        parent_1_str   = "SC5314"
        parent_2_str   = "529L_"
        output_name    = "SCx529L/SCx529L_4way.csv"

    
        

    geno_file = open( geno_file_name, "r" ) # open ze files
    cat_tags_file = open( tags_file_name, "r" )

    geno_file_list = file_reader( geno_file, "\t" ) # Read ze files
    cat_tags_list = file_reader( cat_tags_file, "\t" )

    cat_tag_dict = {}
    for lines in cat_tags_list: # Parse this file list into a dictionary
        if "#" in lines[0]:
            continue
        new_key = lines[2]
        cat_tag_dict[new_key] = lines[3:5]

    geno_file_dict = {}
    for lines in geno_file_list: # Parse this file list into a dictionary too
        new_key = lines[0]
        geno_file_dict[ new_key ] = lines[1:]

    combined_dict = {}
    for cat_nums, mkrs in geno_file_dict.items(): # Combine the two dictionaries

        # I'm gonna take the opposite approach: combine the two dictionaries into a single one and then take what I want out of that; delete nothing! 
        if cat_nums == "# Catalog ID":
            combined_dict[ "# Catalog ID" ] = ["# Catalog ID", "Chromosome", "Chr_Position"] + geno_file_dict[ cat_nums ][3:]
            continue

        combined_dict[ cat_nums ] = [ cat_nums, cat_tag_dict[ cat_nums ][0], cat_tag_dict[ cat_nums ][1] ] + geno_file_dict[ cat_nums ][3:]

    output_file = open( output_name, "w" )
    output_dict = parent_combiner( combined_dict, par_1_str = parent_1_str, par_2_str = parent_2_str )

    printable_list = []
    for mkrs in output_dict.values(): 
        printable_list.append( [ mkrs[1], mkrs[2], mkrs[0] ] + mkrs[3:])

    trp_print_list = transpose( printable_list )
    csv_printer( trp_print_list, output_file )

main(shoot_trouble = 1)
# shoot_trouble = 1 for P6xSC data
# shoot_trouble = 2 for 529LxSC data. 
