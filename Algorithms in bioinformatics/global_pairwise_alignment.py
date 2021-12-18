#!/usr/bin/env python3

"""
Author: Thierry Haddad

Description: This script will perform global, pairwise alignments for given 
sequences. 
"""
#import statements here

# functions between here and __main__
default_penalty = -4 # Default penalty for a gap, overridable.
blosum = """
# http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
      A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
   A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
   R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
   N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
   D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
   C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
   Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
   E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
   H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
   I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
   L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
   K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
   M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
   F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
   P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
   S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
   T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
   W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
   Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
   V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
   B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
   Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
   * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
"""

def blosum62():
    """Return order and similarity scores from BLOSUM62 matrix

    order: dict of {res: idx_in_matrix}
    blosum_matrix: list of lists with similarity scores
    """
    order = {}
    blosum_matrix = []
    for line in blosum.split('\n'):
        if line.startswith('#'):
            continue
        if not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) == 24:
            for idx, sym in enumerate(parts):
                order[sym] = idx
        else:
            # list around the map construction for python3 compatibility
            blosum_matrix.append(list(map(int,parts[1:])))
    return order, blosum_matrix

BLOSUM62_ORDER, BLOSUM62_MATRIX = blosum62()

def score(res1, res2):
    """Return similarity score from BLOSUM62 matrix for two residues
    
    res1: string, amino acid
    res2: string, amino acid
    """
    lookup1 = BLOSUM62_ORDER[res1]
    lookup2 = BLOSUM62_ORDER[res2]
    return BLOSUM62_MATRIX[lookup1][lookup2]

# write your own functions below here

def gap_penalty(gap_length: int) -> int:
    """ Returns gap penatly based on the gap length.
    Applies default gap penatly times the gap length.

    gap_length: int, length of gap in alignment
    """
    penalty = gap_length*default_penalty
    return penalty

def initial_matrix(seq1: str, seq2: str) -> list:
    """ Returns an initial matrix based on two sequences.
    Number of columns is the length of sequence 2, plus 1.
    Number of rows if the length of sequence 1, plus 1.
    Fills the first column and first row with only the gap penalties in dictio-
    ary format.
    The other cells get None.

    matrix: list that contains dictionaries or None, first row/column has 
    gap penalty, rest None.
    """
    columns = len(seq2)+1 # Sequence on top.
    rows = len(seq1)+1 # Sequence on the side.
    # List comprehension that fills the first row and columns with dictionaries
    # that contain gap penalties. Others cells get None.
    matrix = [[{'u':j*default_penalty} if i == 0 else {'l':i*default_penalty} if j == 0 
                    else None for i in range(columns)] for j in range(rows)]
    matrix[0][0] = {'x': 0} # Set starting point to 0.
    return matrix

def global_alignment(matrix: list, seq1: str, seq2: str) -> list:
    for i in range(1, len(matrix)):
        for j in range(1, len(matrix[i])):
            res1 = seq1[i - 1]
            res2 = seq2[j - 1]
            residue_score = score(res1, res2)
            position_score = residue_score

            left = list(matrix[i][j - 1].values())[0] - 4
            up = list(matrix[i - 1][j].values())[0] - 4               
            diagonal = list(matrix[i - 1][j - 1].values())[0] + position_score

            position_score_dict = {'u':up , 'l':left , 'd':diagonal}
            highest_scoring_direction = max(position_score_dict, key=position_score_dict.get)
            position_score = {}
            position_score[highest_scoring_direction] = position_score_dict[highest_scoring_direction]
            matrix[i][j] = position_score
    return matrix

def traceback_matrix(matrix: list, seq1: str, seq2:str) -> str:
    print("\nCalculating traceback...")
    rows = len(matrix)-1
    columns = len(matrix[0])-1
    end_score = 0
    ali_seq1 = ""
    ali_seq2 = ""

    while rows >= 0 and columns >= 0:
        last_position = matrix[rows][columns]
        direction = list(last_position.keys())[0]
        value = last_position[direction]
        end_score += value
        if direction == 'l':
            columns -= 1
            ali_seq1 += " "
            ali_seq2 += seq2[columns]
        elif direction == 'u':
            rows -= 1
            ali_seq2 += " "
            ali_seq1 += seq1[rows]
        elif direction == 'd':
            ali_seq1 += seq1[rows - 1]
            ali_seq2 += seq2[columns - 1]
            rows -= 1
            columns -= 1
        elif direction == 'x':
            rows = -1
            columns = -1
    print("Maximum alignment score in matrix: ", end_score)
    return ali_seq1, ali_seq2

def print_alignment(ali_seq1: str, ali_seq2: str):
    # Reverse sequences to normal direction, due to matrix traceback.
    ali_seq1 = ali_seq1[::-1]
    ali_seq2 = ali_seq2[::-1]

    # Make comparison line between the two sequences.
    # | for a match, * for a mismatch and - for a gap.
    align_rule = ""
    for i in range(0, len(ali_seq1)):
        if ali_seq1[i] == ali_seq2[i]:
            align_rule += "|"
        elif ali_seq1[i] == " " or ali_seq2[i] == " ":
            align_rule += "-"
        else:
            align_rule += "*"
    # Split alignments in 80-char lengths to fit to the screen
    length = 80
    seq1 = [ali_seq1[i:i+length] for i in range(0, len(ali_seq1), length)]
    seq2 = [ali_seq2[i:i+length] for i in range(0, len(ali_seq2), length)]
    ali_line = [align_rule[i:i+length] for i in range(0, len(align_rule), length)]

    # Print alignments per 80-char strings
    print("\nPairwise global alignment:\n")
    print("Legend: | match, * mismatch, - gap\n")
    for i in range(len(seq1)):
        print(seq1[i])
        print(ali_line[i])
        print(seq2[i])
        print("\n")

def print_matrix(matrix: list, seq1: str, seq2: str):
    print("sequence 1: ", seq1)
    print("sequence 2: ", seq2)
    header = "\t"
    for i in seq2:
        header = header + "\t" + i
    print(header)

    seq1 = " " + seq1
    for i in range(len(matrix)):
        row = seq1[i]
        for j in range(len(matrix[i])):
            row = row + "\t" + str(matrix[i][j])
        print(row)

if __name__ == "__main__":
    #default_penalty = -4 # Override default penalty

    seq1 = ("THISLINE")
    seq2 = ("ISALIGNED")

    matrix = initial_matrix(seq1, seq2)
    print("\nInitial matrix:\n")
    print_matrix(matrix, seq1, seq2)

    matrix = global_alignment(matrix, seq1, seq2)
    print("\nGlobal alignment matrix:\n")
    print_matrix(matrix, seq1, seq2)
    ali_seq1, ali_seq2 = traceback_matrix(matrix, seq1, seq2)
    print_alignment(ali_seq1, ali_seq2)

    print("-"*80)
    # seq3: GPA1_ARATH
    seq3 = (
            "MGLLCSRSRHHTEDTDENTQAAEIERRIEQEAKAEKHIRKLLLLGAGESGKSTIFKQIKLLFQTGF"
            "DEGELKSYVPVIHANVYQTIKLLHDGTKEFAQNETDSAKYMLSSESIAIGEKLSEIGGRLDYPRLT"
            "KDIAEGIETLWKDPAIQETCARGNELQVPDCTKYLMENLKRLSDINYIPTKEDVLYARVRTTGVVE"
            "IQFSPVGENKKSGEVYRLFDVGGQRNERRKWIHLFEGVTAVIFCAAISEYDQTLFEDEQKNRMMET"
            "KELFDWVLKQPCFEKTSFMLFLNKFDIFEKKVLDVPLNVCEWFRDYQPVSSGKQEIEHAYEFVKKK"
            "FEELYYQNTAPDRVDRVFKIYRTTALDQKLVKKTFKLVDETLRRRNLLEA"
            )

    # seq4: GPA1_ORYSI
    seq4 = (
            "MGSSCSRSHSLSEAETTKNAKSADIDRRILQETKAEQHIHKLLLLGAGESGKSTIFKQIKLLFQTG"
            "FDEAELRSYTSVIHANVYQTIKILYEGAKELSQVESDSSKYVISPDNQEIGEKLSDIDGRLDYPLL"
            "NKELVLDVKRLWQDPAIQETYLRGSILQLPDCAQYFMENLDRLAEAGYVPTKEDVLYARVRTNGVV"
            "QIQFSPVGENKRGGEVYRLYDVGGQRNERRKWIHLFEGVNAVIFCAAISEYDQMLFEDETKNRMME"
            "TKELFDWVLKQRCFEKTSFILFLNKFDIFEKKIQKVPLSVCEWFKDYQPIAPGKQEVEHAYEFVKK"
            "KFEELYFQSSKPDRVDRVFKIYRTTALDQKLVKKTFKLIDESMRRSREGT"
            )
    matrix = initial_matrix(seq3, seq4)
    matrix = global_alignment(matrix, seq3, seq4)
    ali_seq1, ali_seq2 = traceback_matrix(matrix, seq3, seq4)
    print_alignment(ali_seq1, ali_seq2)
