######################################################################################
#  File            : Pairwise_hw3.py
#  Purpose         : 1> Create GLOBAL alignments between two sequences 
#			using  a dynamic programming alignment algorithm 
#                    2> Create a distance matrix for a set of sequences using 
#	                GLOBAL alignment.
#  Developer       : Han Ngo, Rasmi Lamichhane 
#                    CSCI 4314 Spring 2017, HW3
##############################################################################
#   
#   Sample command line arguments to run the program: 
#   python Pairwise_hw3.py -f large_sample.fasta
#   -f specifies a fasta .txt input file to be read in and analyzed
#
##############################################################################
#
#   Runtime and Space Complexity:
#   	1> To make alignment for two sequences, I use matrix imported from Numpy
#          library, which stores m x n space, given m and n being the length 
#	   of the 2 sequences. Time to memoize and construct all values for all 
#          cells in the table take O(mXn) (as there are two loops being called). 
#	   Back-tracking and constructing the alignment takes O(n),given n 
#	   being the length of the longer sequence in case there are two sequences
#	   of different lengths. 
#
#       2> For comparisons of every single sequence with one another, given k  
#	   being the number of sequences in the list, k*(k-1) tables are made 
# 	   to get all the alignments. Each time a table is completed, use the difference 
#	   score to fill in the cell of the distance matrix. 
#	   O(mXnXkX(k-1)) would be the running time to implement part 2
#	   
##############################################################################

import sys, getopt, math, os.path
import numpy as np
from numpy import array

##############################################################################
# 
#   Specification of the command line to run the program
#               
##############################################################################

def usage():
  print "An input FASTA file (-f)"

##############################################################################
# 
#   Main function to take in argument and manipulate them accordantly
#   When a file is specified, sequences in the file content is read into a list
#   which will be used to passed as parameters to solve two parts of the 
#   homework.
#               
##############################################################################
  
def main(argv):    
    
    try:
        opts, args = getopt.getopt(argv,"hf:")
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
		
	    	# Part 1: Print aligned sequences and the length of the alignment.  
		part1(sequences)
	    	# Part 2: Create alignment matrix and distance matrix		
		part2(sequences)
		
            # If the input file given does not exist, print an error message and exit the program
            else:
                print "ERROR! Input file must exist. Correct usage:"
                usage()
                sys.exit(2)
        
##############################################################################
# 
#   Store all sequences in a list
#   *Sequence names are neglected. Note the case that there are more than 
#   two lines for sequence strings for a sequence name. 
#               
##############################################################################

def readInput(inFile):
    f = open(inFile, 'r')
    lines =  f.read().splitlines()
    temp_seq = ""
    sequences = []    

    for i in range(len(lines)):	
	if lines[i].startswith(">"):
		temp_seq = ""		#renew the temporary sequence to store
	else:
		temp_seq = temp_seq + lines[i]
		try:	
			# check that if new sequence starts the next line, append the last one	
			if lines[i+1].startswith(">"):
				sequences.append(temp_seq)
		except:
			sequences.append(temp_seq)
			
    return sequences

##############################################################################
# 
#   Part 1:   
#   Pair up sequences (sequence 0 with 1, 2 with 3,etc). Neglect the final odd one.
#   Print aligned sequences and the length of the alignment for every two.
#               
##############################################################################

def part1(sequences):
    i = 0
    while i < len(sequences)-1:
	DP_table = makeArray(sequences[i], sequences[i+1]) 
	backTrack(DP_table, sequences[i], sequences[i+1], True)
	i += 2
    return

##############################################################################
# 
#   Part 2: Create alignment matrix and distance matrix
#	Create an alignment matrix (1 vs 2, 1 vs 3, etc, 2 vs 3, etc, N vs N) 
#	using a glocal alignment using the scoring matrix described below.
#	
#	Print distance matrix (all values are between 0-1, use 5 decimal places)
#	meaning (1 - similarity = 1 - #matches/lenght of aligned sequence)
#               
##############################################################################

def part2(sequences): 
    score_matrix = np.zeros((len(sequences), len(sequences)))

    #Set 5 decimals for the division result
    np.set_printoptions(formatter={'float': '{: 0.5f}'.format})

    for i in range(len(sequences)):
	for j in range(len(sequences)):
		if i == j:	#comparing itself, the distance scores 0 
			score_matrix[i][j] = 0
		else:		#comparing with every other sequence, score 1 - similarity as said in header
			DP_table = makeArray(sequences[i], sequences[j])
			score_matrix[i][j] = 1- backTrack(DP_table,sequences[i] ,sequences[j], False)
			    
    print(score_matrix)
    return

##############################################################################
# 
#	Applying Dynamic Programming and Needleman_Wunsh algorithm to 
#       fill in all table cells. 
#               
##############################################################################

def makeArray(seq1, seq2):    
    D = np.zeros((len(seq1)+1, len(seq2)+1), dtype=np.int16)
    
    #Intiation 
    D[1:,0] = range(1,len(seq1)+1)	#first column
    D[0,1:] = range(1,len(seq2)+1)	#first row
    
    #calculate the value for each cells based on values from the top, left, and diagonal
    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            D[i,j] = min(D[i-1, j-1] + getValue(seq1[i-1], seq2[j-1]),
                        D[i-1, j] + 1,
                        D[i, j-1] + 1)
    
    #print(D)
    return D


##############################################################################
# 
#	get subst_score point as 0 as whether the two characters are the same, 
#	1 as they are different. This is to calculate the value from the diagonal  
#               
##############################################################################

def getValue(value_i, value_j): 
    return 0 if value_i == value_j else 1


##############################################################################
# 
#	After finishing the DP table, backtracking is to realign the 
#	two sequences. From the final bottom cell, go back where the value 
#	comes from:
#	Case 1: from diagonal: align 2 characters as match/mismatch
#	Case 2: from left:     gap the first sequence
#	Case 3: from top:      gap the second sequence
#	To benefit part 2 of the homework, this method returns 
#	alignment scores
#               
##############################################################################

def backTrack(D,seq1, seq2, printAlignment):
    aligned_seq1 = ""
    aligned_seq2 = ""
    editDist = ""	#To print out the sequence length after realignment	
    i = len(seq1)
    j = len(seq2)
    score = 0		#To create the distance matrix in part 2	
    
    
    while i>0 and j>0:
	#Priority: check with the diagonal first
        if D[i-1, j-1] + getValue(seq1[i-1], seq2[j-1]) == D[i,j]:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
	  
	    #For part 2: similarity scores 1 	
            if getValue(seq1[i-1], seq2[j-1]) == 0:
	    	score += 1
	    i -= 1
            j -= 1
	
	#Check with the top, plus 1 as the indel cost    
        elif D[i-1, j] + 1 == D[i,j]:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "_" + aligned_seq2
            i -= 1
        
	#Check with the left, 1 is also the indel cost
        else:
            aligned_seq1 = "_" + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
        editDist = "=" + editDist
    
    #In case two sequences are not the same length, finish up the rest
    while i > 0:
            aligned_seq1 = seq1[i-1] + aligned_seq1
	    aligned_seq2 = "_" + aligned_seq2
            editDist = "=" + editDist
            i -= 1
    while j > 0:
	    aligned_seq2 = seq2[j-1] + aligned_seq2
	    aligned_seq1 = "_" + aligned_seq1
            editDist = "=" + editDist
            j -= 1
    
    if printAlignment:		#Just print this out for part 1 only
    	print(aligned_seq1)
    	print(aligned_seq2)
    	print(editDist + " (" + str(len(aligned_seq1)) + ")")

    return score/float(len(aligned_seq1))
    

##############################################################################
# 
#   Trigger the program
#               
##############################################################################

if __name__ == "__main__":
    main(sys.argv[1:])