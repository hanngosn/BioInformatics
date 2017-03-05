######################################################################################
#  File            : HanNgo_hw2.py
#  Purpose         : 1> To find the most common k-mer within the sequences given 
#                    2> To find the k-mer that occurs in the most sequences given 
#                    3> Generate output with GFF format
#  Developer       : Han Ngo 
#                    CSCI 4314 Spring 2017, HW1
##############################################################################
#   
#   Sample command line arguments to run program: 
#   python HanNgo_hw2.py -f sample.fasta.txt -l 4
#   -f specifies a fasta .txt input file to be read in and analyzed
#   -l specifies k-mer length 
#
##############################################################################
#
#   Runtime and Space Complexity:
#   
#   1> Storing file content uses dictionary with keys keeping sequence names 
#      and values keeping sequence strings. Given a FASTA file has n sequences,
#      the runtime takes O(n) for hashing and O(1) for querying sequence. 
#      Space takes up more than n slots according to the default hashing function,
#      but was generally adequately spread out and should not be a big problem with 
#      short-sized DNA and protein.
#   2> Counting sequences and occurances are operated at the same time, taking 
#      O(n*(sequence_length-kmer_length+1)) for the nested loops. Inside the inner loop, 
#      a dictionary will be constructed with keys keeping k-mers and values keeping 
#      occurances and seq counts. Querying information later on will take O(1) also.
#      Space also takes up more than enough kmers with length k. 
#      Benefit: The necessary info is stored in the nested- list dictionary
#   3> Custom sort utilizes the default sorting function in Python with runtime 
#      O(klogk) both on average and worst case, with k being the total number 
#      of kmers with length k.
#
#   Memory limitations: 
#      There are at most 4^k slots with k from 3-8, so memory would consume at 
#      most k*4^k bytes in total (for each character takes up a byte). The 
#      memory limitation is that the slots in dictionary are not fully  
#      occupied because of the default hashing function, so it actually 
#      creates blank space in between, but it should not be a big problem, 
#      cause the querying time with O(1) is more advantageous to us.  
#
##############################################################################

import sys, getopt, math, os.path

##############################################################################
# 
#   Specification of the command line to run the program
#               
##############################################################################

def usage():
  print "An input FASTA file (-f) and size of k-mer (-l) must be specified."

##############################################################################
# 
#   Main function to take in arguments and manipulate them accordantly
#   When a file is specified, the file content is next read into the program memory 
#   through hashing, the specified kmer length is then checked and used as 
#   as a parameter to initiating a dict of kmers with that length. All functions 
#   are executed after these callings and generate output with GFF format.
#               
##############################################################################
  
def main(argv):    
    
    try:
        opts, args = getopt.getopt(argv,"hf:l:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        if opt == '-f':
            if os.path.isfile(arg):
                print "\nFile: " + arg
                sequences = readInput(arg)
            # If the input file given does not exist, print an error message and exit the program
            else:
                print "ERROR! Input file must exist. Correct usage:"
                usage()
                sys.exit(2)
        if opt == '-l':
            try:
                userInput = int(arg)
		
		#Length has to be between 3-8 inclusively
		if userInput < 3 or userInput > 8:
                	print "ERROR! Input number must be from 3 to 8."
                	sys.exit(2)
		print "K-mer length: " + arg
                
		#Start calling functions to manipulate information                
		table = seqCount(sequences,arg)
                list = listConvert(table)
                sorted_occurances = customSort(list, 1)
                sorted_seqCount = customSort(list, 2)
                
                print "\nk-mers\tOccurances"
                output(sorted_occurances, 5, 1)
                print "\nk-mers\tSeq Count"
                output(sorted_seqCount, 5, 2)
            
            #If the input user is not a numerical value, print error message and exit.
            except ValueError:
                print "ERROR! K-mer length must be a digit. Correct usage:"
                usage()
                sys.exit(2)
        
##############################################################################
# 
#   Hashing file content into keys (sequence name) and values (sequence) so that 
#   every time querying takes O(1) 
#               
##############################################################################

def readInput(inFile):
    f = open(inFile, 'r')
    sequenceName = ''
    sequence = ''
    output = {}    

    for line in f:
        line = line.rstrip()    
                
        if line.startswith(">"):        #sequence name starts with ">" in FASTA file
            sequenceName = line
            sequence = ''
        else:       #in case the next-line string is of the same previous sequence 
            sequence = sequence + line
            output[sequenceName[1:]] = sequence
            
    return output

##############################################################################
#     
#    This function seqCount receives 2 parameters, the file content(dict) and
#    kmer length k, and returns a TABLE which keeps the K-MER as KEYS and
#    OCCURANCES & SEQ_COUNT as values. 
#    Each time a kmer is cut out from the sequence, it will be checked if
#    it has already existed in the TABLE or not. 
#    Notice: there is a flag isCounted I used below in order to make sure
#    no matter how many kmers appear in a sequence, the SEQ_COUNT will 
#    be incremented once.
#    
#    This table will have the format like this:
#                {kmerName:[isCounted, occurances, seqCount]}
#               
##############################################################################

def seqCount(content, k):
    
    table = {}    
    
    for eachSeq in content.values():
        for index in range(0,len(eachSeq)-int(k)+1):
            isCounted = False                        
            kmer = eachSeq[index:index+int(k)]       #cut out k-mer with length k

            
            if kmer in table:           #check whether the kmer already exists in the table
                table[kmer][1] += 1     #increment occurances
                if not table[kmer][0]:
                    table[kmer][2] += 1             
                    table[kmer][0] = True
        
            else:
                table.setdefault(kmer, [])      
                table[kmer].append(True)            
                table[kmer].append(1)               #first-time occurance
                table[kmer].append(1)               #first-time seqCount
    
        for eachkmer in table:          #end of each sequence, reset to false meaning 
            table[eachkmer][0] = False  #the seqCount has not been incremented yet
                    

#     print table
    return table

##############################################################################
# 
#    Sorting the list within the values of the dictionary gets hard, so 
#    this function listConvert will turn the dictionary into a nested list with
#    each element keeping k-mer name, occurances, and seqCount. 
#               
##############################################################################

def listConvert(dic):
    list = []
    for each in dic.keys():
        sublist = [each, dic[each][1], dic[each][2]]
        list.append(sublist)
        
#    print(list)
    return list

##############################################################################
# 
#    This function uses Python's default bubble sort with my customization.
#    It will be called two times, one to sort the occurances and the other to
#    sort the seqCount.
#
##############################################################################

def customSort(list,column):
    
    return sorted(list, key=lambda x: int(x[column]), reverse = True)

##############################################################################
# 
#   The requirement is to print top 5 kmer with their counts, but if the 
#   fifth has several common ones, print them all out.  
#
##############################################################################

def output(list, limit, column):
    
    index = 0
     
    while index < limit:    
        print(list[index][0] + "\t" + str(list[index][column]))
        try:
	    if index == limit - 1: 
            	if list[index][column] == list[index+1][column]: #compare if their counts are unique 
			limit += 1
        except: 
            break
	
	index += 1
    
    return
    
##############################################################################
# 
#   Program-triggering 
#               
##############################################################################

if __name__ == "__main__":
    main(sys.argv[1:])