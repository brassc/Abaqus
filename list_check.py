# Read the contents of the two text files and split them into lists of numbers.
# List 1 is the big list, List 2 is the smaller list inside list 1. This program creates
# List1-list2=newlist 
# e.g. list_check(list1='GM_Node_Set_NodeList.txt', list2='Swelling_Region_DC_Side_NodeList.txt', filename_list='node_numbers.txt')

def list_check(list1, list2, filename_list):
    with open(list1, 'r') as file1:
        list1 = [int(number) for number in file1.read().split(',')]

    with open(list2, 'r') as file2:
        list2 = [int(number) for number in file2.read().split(',')]

    # Find numbers from list1 that are not in list2
    not_in_list2 = [number for number in list1 if number not in list2]

    # Save the numbers not in list2 to a new text file
    with open(filename_list, 'w') as output_file:
        output_file.write(','.join(map(str, not_in_list2)))

    print("Numbers not in list2 have been saved to '{}'.".format(filename_list))