from Bio import SeqIO
from models import linkedNode


def createMatrix(sequence1, sequence2):
    matrix = []
    for i in range(len(sequence1) + 1):
        matrix.append([])
    for i in range(len(sequence1) + 1):
        for j in range(len(sequence2) + 1):
            matrix[i].append(linkedNode(0, 0, "NULL"))

    return matrix


def openFiles(user_input, success):
    # get user input for fasta files
    while success is False:
        # if user_input is 1, then only open one file and read the first two sequences from that file
        if user_input == "1":
            # make sure the user is inputing the appropriate file types
            try:
                file1 = input("Name or path of FASTA file: ")
                if ".fasta" not in file1:
                    raise OSError("[Custom Error] Incorrect file type.")
                with open(file1) as fasta_file:
                    sequences = SeqIO.parse(fasta_file, 'fasta')
                    sequences = list(sequences)

                sequence1 = str(sequences[0].seq)
                # print(sequence1)
                sequence2 = str(sequences[1].seq)
                # print(sequence2)
                success = True
            except OSError as e:
                print(e)
            except IndexError as ie:
                print(ie)
                print("You provided a file with only one sequence instead of two")

        # if user_input is 2, then open both files and read the first sequence from both files
        elif user_input == "2":
            # make sure the user is inputing the appropriate file types
            try:
                file1 = input("Name or path of FASTA file #1: ")
                if ".fasta" not in file1:
                    raise OSError("[Custom Error] Incorrect file type.")
                with open(file1) as fasta_file:
                    sequences1 = SeqIO.parse(fasta_file, 'fasta')
                    sequences1 = list(sequences1)
                sequence1 = str(sequences1[0].seq)
                # print(sequence1)

                file2 = input("Name or path of FASTA file #2: ")
                if ".fasta" not in file2:
                    raise OSError("[Custom Error] Incorrect file type.")
                with open(file2) as fasta_file:
                    sequences2 = SeqIO.parse(fasta_file, 'fasta')
                    sequences2 = list(sequences2)
                sequence2 = str(sequences2[0].seq)
                # print(sequence2)
                success = True
            except OSError as e:
                print(e)
            except IndexError as ie:
                print(ie)
                print("You provided a file with no sequences")

    return [sequence1, sequence2]


def findOptimalScore(left, above, diagonal, seq1, seq2):
    """
    given three nodes (above, diagonal, and left), this function will find the optimal score and return
    that value to the current node.
    """
    matchMismatch = 0

    if seq1 == seq2:
        matchMismatch = 5  # the nucleotides were a match
    else:
        matchMismatch = -1  # the nucleotides were a mismatch

    option1 = diagonal + matchMismatch  # add diagonal value to mismatch value
    option2 = above - 4  # gap penalty
    option3 = left - 4  # gap penalty

    # find the max of the three options
    result = max(option1, option2, option3)

    # create a new 'resultNode' based of what the result is
    if result == option1:
        resultNode = linkedNode(result, diagonal, "diagonal")
    elif result == option2:
        resultNode = linkedNode(result, above, "above")
    else:
        resultNode = linkedNode(result, left, "left")

    return resultNode


def globalAlignmentMatrix(matrix, seq1, seq2):
    """
    given an empty matrix and two sequences, this function will fill out the matrix with the optimal scores
    """
    sourceNode = linkedNode(0, 0, "NULL")  # source node is always zero
    matrix[0][0] = sourceNode  # initialize source node

    # fill out the matrix
    for i in range(len(seq1) + 1):
        for j in range(len(seq2) + 1):

            # fill the first row with gap penalties
            if i is 0 and j > 0:
                newNode = linkedNode(matrix[i][j - 1].value - 4, matrix[i][j - 1].value, "left")
                matrix[i][j] = newNode
            # fill the first column with gap penalties
            elif j is 0 and i > 0:
                newNode = linkedNode(matrix[i - 1][j].value - 4, matrix[i - 1][j].value, "above")
                matrix[i][j] = newNode
            # fill the rest of the spots in the matrix with the optimum scores
            elif i > 0 and j > 0:
                newNode = findOptimalScore(matrix[i][j - 1].value, matrix[i - 1][j].value, matrix[i - 1][j - 1].value,
                                           seq1[i - 1], seq2[j - 1])
                matrix[i][j] = newNode

    return matrix


def getOptimalAlignment(matrix, seq1, seq2):
    """
    given a matrix and two sequences, this function will return three strings
    line1: the optimal alignment to the second sequence for the first sequence
    line2: shows if there is a match, mismatch, or gap between line1 and line3
    line3: the optimal alignment to the first sequence for the second sequence

    this function backtracks from the very last value in the matrix
    """
    line1 = ""
    line2 = ""
    line3 = ""
    index1 = len(seq1)
    index2 = len(seq2)
    currentNode = matrix[len(seq1)][len(seq2)]

    while index1 != 0 or index2 != 0:
        # if the current nodes source direction is 'diagonal', there was a match or mismatch
        if currentNode.direction == "diagonal" and seq1[index1 - 1] == seq2[index2 - 1]:  # match
            line1 = line1 + seq1[index1 - 1]
            index1 -= 1
            line2 = line2 + "*"
            line3 = line3 + seq2[index2 - 1]
            index2 -= 1
            currentNode = matrix[index1][index2]
        # if the current nodes source direction is 'diagonal', there was a match or mismatch
        elif currentNode.direction == "diagonal" and seq1[index1 - 1] != seq2[index2 - 1]:  # mismatch
            line1 = line1 + seq1[index1 - 1]
            index1 -= 1
            line2 = line2 + "|"
            line3 = line3 + seq2[index2 - 1]
            index2 -= 1
            currentNode = matrix[index1][index2]
        # if the source nodes direction is 'left', put a gap in the first sequence
        elif currentNode.direction == "left":
            line1 = line1 + "-"
            line2 = line2 + " "
            line3 = line3 + seq2[index2 - 1]
            index2 -= 1
            currentNode = matrix[index1][index2]
        # if the source nodes direction is 'above', put a gap in the second sequence
        elif currentNode.direction == "above":
            line1 = line1 + seq1[index1 - 1]
            index1 -= 1
            line2 = line2 + " "
            line3 = line3 + "-"
            currentNode = matrix[index1][index2]

    return [line1, line2, line3]
