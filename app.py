
import time
from functions import getOptimalAlignment, globalAlignmentMatrix, openFiles, createMatrix

user_input = ""
success = False


if __name__ == '__main__':

    # get user input
    # if user input is invalid, loop will continue
    while user_input != "1" and user_input != "2":
        user_input = input("Number of FASTA files: ")
        if user_input != "1" and user_input != "2":
            print("invalid input")

    # get sequences
    sequences = openFiles(user_input, success)
    sequence1 = sequences[0]
    sequence2 = sequences[1]

    # convert sequences to upper case
    sequence1 = sequence1.upper()
    sequence2 = sequence2.upper()

    # print sequences
    print("")
    print("Sequence 1: " + sequence1)
    print("Sequence 2: " + sequence2)
    print("")
    time.sleep(3)

    # create matrix of null linkedNode objects with dimensions m x n
    # where m is len(sequence1) + 1 and n is len(sequence2) + 1
    matrix = createMatrix(sequence1, sequence2)

    # get filled out matrix from globalAlignmentMatrix() function defined in functions.py
    matrix = globalAlignmentMatrix(matrix, sequence1, sequence2)

    # get optimal alignment from getOptimalAlignment() function defined above
    opt_alignment = getOptimalAlignment(matrix, sequence1, sequence2)

    # output the match bonus, mismatch penalty, and gap penalty
    print("Match: [+5]     Mismatch: [-1]     Gap Penalty: [-4]")
    print("")
    time.sleep(3)

    # print the filled out matrix with formatting
    print("Optimal Global Alignment Matrix:")
    print("")
    for row in matrix:
        print("[", end=" ")
        for val in row:

            if len(str(val.value)) == 1:
                print(" " + str(val.value), end="  ")
            elif len(str(val.value)) == 2:
                print(str(val.value), end="  ")
            else:
                print(str(val.value), end=" ")
        print("]\n")
        time.sleep(1)

    # optimal global alignment score is the very last value of the matrix (bottom right)
    print("The optimal global alignment score is " + str(matrix[len(sequence1)][len(sequence2)].value))
    print("")
    time.sleep(3)

    # print the optimal global alignment
    print("The optimal global alignment is: ")
    print("")
    print(opt_alignment[0][::-1])  # print backwards
    print(opt_alignment[1][::-1])
    print(opt_alignment[2][::-1])
    print("\n")
    print("[ - ]     GAP")
    print("[ * ]     MATCH")
    print("[ | ]     MISMATCH")

