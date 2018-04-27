#!/usr/bin/env python
import sys
import string
#open file                                                                                      <<<<<<<<<<<<<<<<<<<<<<<<<
filename = "3fasta.fa" # you can change the input/path/file here!       <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
C=seq[2].strip()
#A="ATGGCCG"                                        # you can uncomment these
#B="AGTAAATTG"                                    # three lines to be independent from file input
#C="ATGTAAC"                                         #
gap = -2
match = 1
mismatch = -1
def scoring(a,b):
    if a==b:
        score = 2
        return score
    if a!=b:
        score = -2
        return score
def score_3seq(a,b,c):
    if a == b == c: return 2+2
    if a != b == c: return 2-2
    if a == b != c: return 2-2
    if a != b != c: return -2-2

#make similarity matrix
def similarity(A, B, C):
    len_a = len(A) #s #n
    len_b = len(B) #t #m
    len_c = len(C)
    #initialize matrix
    matrix = [[[0 for x in range(len_c+1)] for y in range(len_b+1)]for z in range(len_c+1)]

    #initialize edges
    for i in range(1, len_a + 1):
        matrix[i][0][0] = i * gap
    for j in range(1, len_b + 1):
        matrix[0][j][0] = j * gap
    for k in range(1, len_c + 1):
        matrix[0][0][k] = k * gap

    #initialize planes
    k=0
    for i in range(1,len_a + 1):
        for j in range(1, len_b+1):
            case_1 = matrix[i - 1][j][k] + gap
            case_2 = matrix[i][j - 1][k] + gap
            case_3 = matrix[i - 1][j - 1][k] + scoring(A[i-1], B[j-1])
            matrix[i][j][k] = max(case_1, case_2, case_3)
    i=0
    for j in range(1, len_b + 1):
        for k in range(1, len_c + 1):
            case_1 = matrix[i][j - 1][k] + gap
            case_2 = matrix[i][j][k - 1] + gap
            case_3 = matrix[i][j - 1][k - 1] + scoring(B[j-1], C[k-1])
            matrix[i][j][k] = max(case_1, case_2, case_3)
    j=0
    for i in range(1, len_a + 1):
        for k in range(1, len_c + 1):

            case_1 = matrix[i-1][j][k] + gap
            case_2 = matrix[i][j][k - 1] + gap
            case_3 = matrix[i-1][j][k - 1] + scoring(A[i - 1], C[k - 1])
            matrix[i][j][k] = max(case_1, case_2, case_3)

    #now three dimensions
    exrapoint=1
    for i in range(1, len_a + 1):
        for j in range(1, len_b + 1):
            for k in range(1, len_c + 1):
                case_1 = matrix[i - 1][j][k] + gap
                case_2 = matrix[i][j][k - 1] + gap
                case_3 = matrix[i][j-1][k] + gap
                case_4 = matrix[i - 1][j - 1][k] + scoring(A[i - 1], B[j - 1])
                case_5 = matrix[i - 1][j][k - 1] + scoring(A[i - 1], C[k - 1])
                case_6 = matrix[i][j - 1][k - 1] + scoring(B[j - 1], C[k - 1])
                case_7 = matrix[i - 1][j-1][k - 1] + score_3seq(A[i - 1], B[j - 1], C[i-1])
                matrix[i][j][k]=max(case_1, case_2, case_3, case_4, case_5, case_6, case_7)
    return matrix
M=similarity(A,B,C)
#for plane in M:
#    for line in plane:
#        print(line)

#Backtracking
aaa=[]
bbb=[]
ccc=[]
def align(i,j, k):
    if i==0 and j==0 and k==0:
        pass
    elif M[i][j][k] == M[i-1][j][k] + gap:
        align(i-1, j,k)
        aaa.append(A[i-1])
        bbb.append('-')
        ccc.append('-')
    elif M[i][j][k] == M[i][j - 1][k] + gap:
        align(i, j-1, k)
        aaa.append('-')
        bbb.append(B[j-1])
        ccc.append('-')
    elif M[i][j][k] == M[i][j][k - 1] + gap:
        align(i, j, k-1)
        aaa.append('-')
        bbb.append('-')
        ccc.append(C[j-1])
    elif M[i][j][k] == M[i - 1][j - 1][k] + scoring(A[i - 1], B[j - 1]):
        align(i - 1, j - 1, k)
        aaa.append(A[i - 1])
        bbb.append(B[j - 1])
        ccc.append('-')
    elif M[i][j][k] == M[i][j - 1][k - 1] + scoring(B[j - 1], C[k - 1]):
        align(i, j - 1, k - 1)
        aaa.append('-')
        bbb.append(B[j - 1])
        ccc.append(C[k - 1])
    elif M[i][j][k] == M[i - 1][j][k - 1] + scoring(A[i - 1], C[k - 1]):
        align(i - 1, j, k - 1)
        aaa.append(A[i - 1])
        bbb.append('-')
        ccc.append(C[k - 1])
    elif M[i][j][k] == M[i - 1][j - 1][k - 1] + score_3seq(A[i - 1], B[j - 1], C[i - 1]):
        aaa.append(A[i - 1])
        bbb.append(B[j - 1])
        ccc.append(C[k - 1])
    return(aaa, bbb, ccc)

alignment = align(len(A),len(B), len(C))

print('sequence 1: {}\nSequence 2: {}\nSequence 2: {}\nScore: {}\nAlignment: \n'
      '{}\n{}\n{}'.format(A, B, C, M[len(A)][len(B)][len(C)], alignment[0], alignment[1], alignment[2]))