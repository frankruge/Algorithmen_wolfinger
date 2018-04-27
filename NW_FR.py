#!/usr/bin/env python
import sys
import string
#open file                                                                                      <<<<<<<<<<<<<<<<<<<<<<<<<
filename = "2fasta.fa" # you can change the input/path/file here!       <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
textfile = open(filename, 'r') #                                                                <<<<<<<<<<<<<<<<<<<<<<<<<
message = filename+" evaluated:"
if len(sys.argv)>1: # This program is designed to get file input from the command line
    try:
        textfile = open(sys.argv[1], "r")
        message = sys.argv[1]+" evaluated:"
    except:
        textfile = open(filename, 'r')
        message = sys.argv[1]+" does not seem to exist. used file ./rdm50/rmd50_11 instead :"
fastalist=[]
seq=[]
count=0
for line in textfile:
    if not line.startswith('>'):
        seq.append(line)
A=seq[0].strip()
B=seq[1].strip()
#A="ATGGCCG"                                        # you can uncomment these
#B="AGTAAAGGCCG"                                    # two lines to be independent from file input
gap = -2
match = 1
mismatch = -1

#make similarity matrix
def similarity(A, B):
    len_a = len(A) #s #n
    len_b = len(B) #t #m
    #initialize matrix
    matrix = [[0 for x in range(len_b+1)] for y in range(len_a+1)]
    for i in range(1, len_a+1):
        matrix[i][0] = i * gap
    for i in range(1, len_b+1):
        matrix[0][i] = i * gap
    for i in range(1,len_a+1):
        for j in range(1, len_b+1):
            case_1 = matrix[i-1][j] + gap #sequence A opens gap
            case_2 = matrix[i][j-1] + gap #sequence B opens gap
            if A[i-1] == B[j-1]:
                penalty = match
            elif A[i-1] != B[j-1]:
                penalty = mismatch
            case_3 = matrix[i-1][j-1] + penalty
            matrix[i][j]=max(case_1, case_2, case_3)
    return matrix
M=similarity(A,B)
for line in M:
    print(line)

#Backtracking
aaa=[]
bbb=[]
def align(i,j):
    if i==0 and j==0:
        pass
    elif M[i][j] == M[i-1][j] + gap:
        align(i-1, j)
        aaa.append(A[i-1])
        bbb.append('-')
    elif M[i][j] == M[i][j-1] + gap:
        align(i, j-1)
        aaa.append('-')
        bbb.append(B[j-1])
    elif M[i][j] == M[i-1][j-1] + match or M[i][j] == M[i-1][j-1] + mismatch:
        align(i - 1, j - 1)
        aaa.append(A[i - 1])
        bbb.append(B[j - 1])
    return(aaa, bbb)

alignment = align(len(A),len(B))

print('sequence 1: {}\nSequence 2: {}\nScore: {}\nAlignment: \n'
      '{}\n{}'.format(A, B, M[len(A)][len(B)], alignment[0], alignment[1]))