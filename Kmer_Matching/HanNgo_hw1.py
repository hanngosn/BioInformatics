######################################################################################
#  File            : HanNgo_hw1.py
#  Purpose         : To calculate the GC content of nucleotide sequences and find
#                    exact matches for an input k-mer in nucleotide and protein
#                    sequences	
#  Developer       : Han Ngo 
#		     CSCI 4314 Spring 2017, HW1
##############################################################################
#   
#   Sample command line arguments to run program: 
#   python HanNgo_hw1.py -f sample.fasta.txt -s species001 -k ATG
#   -f specifies a fasta .txt input file to be read in and analyzed
#   -s specifies the sequence to be evaluated
#   -k specifies the k-mer to be compared
#
##############################################################################
#
#   Runtime and Space Complexity:
#
#   Runtime:
#     The algorithm to find GC content run theta(m) with size of specified 
#     sequence given m. The algorithm to find matching positions between k-mer
#     and the specified sequence took theta(m-k+1) to loop through all the possible 
#     character in the sequence with k length.
#   Space:
#     The algorithm utilizes hashmap/dictionary to keep sequence names as keys
#     and sequence strings as values. The space may be taken up more than the 
#     total number of sequence 
#
##############################################################################
# 
#   References: Arguments manipulation part from Alex Okeson HW2 2016
#               
##############################################################################
import sys, getopt, math, os.path

##############################################################################
# 
#   Specification of the command line to run the program
#               
##############################################################################
def usage():
  print "An input file (-f), sequence name (-s), and kmer (-k) must be specified."

##############################################################################
# 
#   Main function to take in arguments and manipulate them accordantly
#   When a file is specified, the file content is next read into the program memory 
#   through hashing, the specifid sequence name is then checked in the hash keys 
#   and gives back the value of the string sequence. K-mer is then manipulated into
#   its reverse complement and matching position is then checked. After all, a file 
#   format GFF is generated.
#               
##############################################################################
  
def main(argv):	
	seqName = ''
	seq = ''
	forward = True 

	try:
		opts, args = getopt.getopt(argv,"hf:s:k:")
	except getopt.GetoptError:
		usage()
		sys.exit(2)
	
	for opt, arg in opts:
		if opt == '-h':
			print "Correct usage is:"
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
		if opt == '-s':
			if sequences.has_key(arg):	
				seqName = arg
				seq = sequences[arg]
				print "Sequence name found: " + seqName + "\n" + "Sequence: " + seq + "\nSequence length: " + str(len(seq)) + "\nGC Content: " + str(calculateGC(seq)) + "%"	 
			else:
				print "ERROR! Sequence name must exist. Correct usage:"
				usage()
				sys.exit(2)
		if opt == '-k':
			print "K-mer: " + arg + "\n"
			try: 	
				id = 0
				id = generateGFF(seqName, seq, kmerMatch(seq, arg), arg, forward,id)
				
				forward = False										
				generateGFF(seqName, seq, kmerMatch(seq, reverseComp(arg)), arg, forward,id)

			except:
				print "No matches found"	
			

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
		if line.startswith(">"):
			sequenceName = line
			sequence = ''
		else:
			sequence = sequence + line
			output[sequenceName[1:]] = sequence
	
	#print output
	return output

##############################################################################
# 
#   Calculate the content of GC by looping through each letter in the sequence 
#               
##############################################################################

def calculateGC(seq):
	g = 0
	a = 0
	c = 0
	t = 0

	for char in seq:
		if char == "G":
			g += 1
		if char == "A":
			a += 1
		if char == "C":
			c += 1
		if char == "T":
			t += 1
	#+0. to get the float value 
	gc = (g+c+0.)/(g+a+c+t+0.)*100
	
	#print gc	
	return gc	

##############################################################################
# 
#   This function accepts a k-mer as a parameter and return the reverse complement 
#   of it
#               
##############################################################################

def reverseComp(kmer):
	reverse = '' 
	for each in kmer:
		if each == 'A':
			reverse = 'T'+reverse
		if each == 'T':
			reverse = 'A'+reverse
		if each == 'G':
			reverse = 'C'+reverse
		if each == 'C':
			reverse = 'G'+reverse		
	
	#print reverse
	return reverse				

##############################################################################
# 
#   This function returns all the positions that the k-mer matches the sequence
#   by looping through all the characters and comparing each k-length subsequence
#   with the k-mer
#               
##############################################################################

def kmerMatch(seq, kmer):
	index = 0
	matchPos = []

	while (index < len(seq)-len(kmer)+1):
		if kmer == seq[index:index+len(kmer)]:
			matchPos.append(index)					
		index += 1	
	
	#print matchPos 
	return matchPos
	

##############################################################################
# 
#   This is to print out to terminal the inquiries under GFF format, in 
#   accordance with its specification like strand +/- for forward/backward, id
#   number for each, starting and ending position of the feature.
#               
##############################################################################

def generateGFF(seqName, seq, posMatch, kmer, forward, id):

	for index in range(0,len(posMatch)):
		start = posMatch[index] + 1
		end = posMatch[index] + len(kmer) 
		if forward:
			print seqName + "\tHannie\tMatch\t" + str(start) + "\t" + str(end) + "\t100.\t+\t.\tID=" + str(id+index+1)  	
		else:
			print seqName + "\tHannie\tMatch\t" + str(start) + "\t" + str(end) + "\t100.\t-\t.\tID=" + str(id+index+1)  	
	return id+index+1		

##############################################################################
# 
#   Program-triggering 
#               
##############################################################################

if __name__ == "__main__":
	main(sys.argv[1:])
