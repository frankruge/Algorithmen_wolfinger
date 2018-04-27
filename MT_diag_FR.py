#!/usr/bin/env python
import sys

#open file                                                                                      <<<<<<<<<<<<<<<<<<<<<<<<<
filename = "./rdm50/rmd50_11" # you can change the input/path/file here!       <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
textfile = open(filename, 'r') #                                                                <<<<<<<<<<<<<<<<<<<<<<<<<
message = filename+" evaluated:"
if len(sys.argv)==2: # This program is designed to get file input from the command line
    try:
        textfile = open(sys.argv[1], "r")
        message = sys.argv[1]+" evaluated:"
    except:
        textfile = open(filename, 'r')
        message = sys.argv[1]+" does not seem to exist. used file ./rdm50/rmd50_11 instead :"

#read file into lists
gright=[]
gdown=[]
gdiag=[]
switch = "down"
for line in textfile:
    if line.startswith('  ') and switch == "down":
        gdown.append(line.strip().split())
        #print(switch)
    if line.startswith("---")and switch == "down":
        switch = "right"
        continue
        #print(switch)
    if line.startswith('  ') and switch == "right":
        gright.append(line.strip().split())
        #print(switch)
    if line.startswith("---")and switch == "right":
        switch = "diag"
        continue
        #print(switch)
    if line.startswith('  ') and switch == "diag":
        gdiag.append(line.strip().split())
        #print(switch)
    if line.startswith('G'):
        continue
textfile.close()

# print the three lists
print("gdown")
for line in gdown:
	print(line)
print("#################################################################")
print("gright")
for line in gright:
	print(line)
print("#################################################################")
print("gdiag")
for line in gdiag:
	print(line)

def manhattan_tourist(down, right, diag):
    #initialize scoring matrix
    cols=len(down[0])
    rows=len(right)
    scoring_matrix = [[0 for x in range(cols)] for y in range(rows)] #initialized with 0 values
    backtrack_m= [[0 for x in range(cols)] for y in range(rows)] #initialization with 0 for debugging

    # fill up first row - scoring_matrix[0][i]
    for i in range(1,cols):
        scoring_matrix[0][i] = round(float(scoring_matrix[0][i - 1]) + float(right[0][i - 1]), 2)
        backtrack_m[0][i] = '-'

    # fill up first column - scoring_matrix[1][0]
    for i in range(1, rows):
        scoring_matrix[i][0] = round(float(scoring_matrix[i-1][0]) + float(down[i-1][0]), 2)
        backtrack_m[i][0] = '-'

    # calculate scores including diagonal values
    for i in range(1,rows):
        for j in range(1,cols):
            scoring_matrix[i][j] = round(max((scoring_matrix[i-1][j] + float(down[i-1][j])),    #  "↓" previous row; this column
                    (scoring_matrix[i][j-1] + float(right[i][j-1])),                            #  "→" this row; previous column
                    (scoring_matrix[i-1][j-1] + float(diag[i-1][j-1]))), 2)                     #  "↘" previous row; previous column
    # backtracking matrix
    for i in range(rows):
        for j in range(cols):
            if scoring_matrix[i][j] == round(scoring_matrix[i-1][j-1] + float(diag[i-1][j-1]), 2):
                backtrack_m[i][j]= "↘"
            if scoring_matrix[i][j] == round(scoring_matrix[i-1][j] + float(down[i-1][j]), 2):
                backtrack_m[i][j]= "↓"
            if scoring_matrix[i][j] == round(scoring_matrix[i][j-1] + float(right[i][j-1]), 2):
                backtrack_m[i][j]= "→"

    #print scoring_matrix and backtrack_m
    print("#################################################################")
    for line in scoring_matrix:
        print(line)
    print("#################################################################")
    for line in backtrack_m:
        print(line)


    return scoring_matrix[rows-1][cols-1]

a = manhattan_tourist(down=gdown, right=gright, diag=gdiag)
print(message)
print("maximum score is "+str(a))
