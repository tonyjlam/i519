'''importing argparse'''
from argparse import ArgumentParser

'''biopython'''
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#open fasta and determine nucleotide freq

def calc_comp(sequence):

	"""calc_comp annalyzes composition of nucleotide composition
           adds composition ratio into dictionary
	   returns dictionary

	Args:
	   'fasta file'
	Returns:
 	   dictionary of each nt sequence"""
	
	#empty dictionary
	comp = {}	
	
	#accesses sequence from fasta input
	sequence1 = sequence.seq

	#dictionary entries for nt count	
	seqlen = len(sequence1)
	comp['A'] = (float(sequence1.count('A'))/seqlen)
	comp['G'] = (float(sequence1.count('G'))/seqlen)
	comp['T'] = (float(sequence1.count('T'))/seqlen)
	comp['C'] = (float(sequence1.count('C'))/seqlen)
	
	return comp

from numpy.random import choice
'''imports choice from numpy.random '''

def gen_sequence(comp, length):

	'''takes nucleotide composition and integer and returns random sequence
	with nt composition and  the integer as the random sequence length

	Args:
	   comp: composition returned from function calc_comp
	   length: desired length of the random sequence
	Returns:
	   returns a random sequence with given composition and length.'''

	#empty storage string
	seq_string = ''

	#possible nucleotides
	nt = ['A', 'C', 'T', 'G']
	
	#assigns generated list to ntlist
	#uses numpy's weighted probability to randomly pick a nt	
	ntlist =  choice(nt, length, p=[comp['A'], comp['C'], comp['T'], comp['G']])

	#joins seq_string and makes list into string
	seq_string = ''.join(ntlist)

	
	sequence = Seq(seq_string)
	return SeqRecord(sequence, id='Random Sequence', description=comp.__repr__())

#ALL MODIFICATIONS GO ABOVE THIS LINE #####

### Part 0 : Argument Parsing
### We want out program to have easy-to-use parameters
### We are using the argparse library for this

parser = ArgumentParser(description='Reads a fasta sequence and generates random ones with the same compositon')
# input file
parser.add_argument('-i', '--infile', help='input file in fasta format')
# output file
parser.add_argument('-o', '--outfile', help='output file in fasta format')
# number of samples (default is 1)
parser.add_argument('-n', '--nrandom', type=int, default=1, help='number of random sequences to be generated')

args=parser.parse_args()

### Part 1 : Reading a FASTA file
### Biopython library makes this extremely easy using the SeqIO module

print 'Parsing the input sequence...'
input_sequence = SeqIO.read(args.infile, 'fasta')
# Note that this returns a SeqRecord object native to Biopython

# we can easily get the length of sequence
n_in_seq = len(input_sequence)

### Part 2 : Getting the nucleotide composition of the input sequence
print 'Calculating input sequence composition...'
input_comp = calc_comp(input_sequence)

### Part 3 : Generating the random sequences with the same composition
print 'Generating random sequences...'

# first define an empty list
random_sequences = list()

# then use a for loop to populate the list
for i in range(args.nrandom) :
    random_sequences.append( gen_sequence(input_comp, n_in_seq) )

# Here is 'the python way' of doing the same thing in a single line
# This uses 'list comprehension' which is a very versatile feature of python
# random_sequences = [ gen_sequence(input_comp, n_in_seq)) for i in range(args.nrandom) ]

### Part 4 : Writing random sequences to a FASTA file
print 'Writing sequences to file...'
SeqIO.write(random_sequences, args.outfile, 'fasta')
