#!/usr/bin/python3

"""
DESCRIPTION:
    Template code for the FIRST Advanced Question of the Hidden Markov Models
    assignment in the Algorithms in Sequence Analysis course at the VU.

INSTRUCTIONS:
    Complete the code (compatible with Python 3!) upload to CodeGrade via
    corresponding Canvas assignment. Note this script will be graded manually,
    if and only if your "hmm.py" script succesfully implements Baum-Welch
    training! Continuous Feedback will not be available for this script.

AUTHOR:
    Janneke Hulsen (2751395)
"""

from argparse import ArgumentParser, RawTextHelpFormatter
from hmm_utility import load_tsv
from numpy.random import choice



def parse_args():
    #####################
    # START CODING HERE #
    #####################
    # Implement a simple argument parser (WITH help documentation!) that parses
    # the information needed by main() from commandline. Take a look at the
    # argparse documentation, the parser in hmm_utility.py or align.py
    # (from the Dynamic Programming exercise) for hints on how to do this.

    "Parse inputs from commandline and return them as Namespace objects."

    parser = ArgumentParser(prog = 'python3 sequence_generator.py',
    formatter_class = RawTextHelpFormatter, description =
    '  Generate sequences from given transmission and emission matrices.\n\n'
    '  Example syntax:\n'
    '    python3 sequence_generator.py A.tsv E.tsv')

    # Positional arguments.
    parser.add_argument('transition', help='path to a TSV formatted transition matrix')
    parser.add_argument('emission', help='path to a TSV formatted emission matrix')

    # Optional arguments.
    parser.add_argument('-n', dest = 'seqnumber', type = int, default = 12,
                       help = 'number of sequences to generate')
    parser.add_argument('-v', '--verbose', dest='verbosity', action='count', default=0,
        help='print verbose output\n')
    parser.add_argument('-o', dest='out_dir', type = str, default = "Generated_sequences.fasta",
        help='path to a directory where output file is saved\n'
             '  (directory will be made if it does not exist)')

    return parser.parse_args()

    #####################
    #  END CODING HERE  #
    #####################


def generate_sequence(A,E):
    #####################
    # START CODING HERE #
    #####################
    # Implement a function that generates a random sequence using the choice()
    # function, given a Transition and Emission matrix.

    allStates = list(A.keys())
    emittingStates = allStates[1:len(allStates)-1] # Remove Begin and End state.

    # Initialize.
    stateseq = "B"

    # Add states based on transition matrix.
    i = 0
    sequence = ""
    while stateseq[-1] != "E":
        # Add next state.
        probabilities = []
        for key in A[stateseq[i]].keys():
            probabilities.append(A[stateseq[i]][key])
        stateseq += list(choice(allStates, 1, p = probabilities)[0])[0]

        # Add next amino acid.
        probabilities = []
        aa = []
        if stateseq[i] in emittingStates:
            for key in E[stateseq[i]].keys():
                aa.append(key)
                probabilities.append(E[stateseq[i]][key])
            sequence += list(choice(aa, 1, p = probabilities)[0])[0]

        # Next position.
        i += 1
        
    #####################
    #  END CODING HERE  #
    #####################

    return sequence



def main():
    args = parse_args()
    #####################
    # START CODING HERE #
    #####################
    # Uncomment and complete (i.e. replace '?' in) the lines below:

    N = args.seqnumber               # The number of sequences to generate
    out_file = args.out_dir        # The file path to which to save the sequences
    A = load_tsv(args.transition)    # Transition matrix
    E = load_tsv(args.emission)    # Emission matrix
    verbosity = args.verbosity    # To print or not to print

    with open(out_file, "w") as f:
        for i in range(N):
            seq = generate_sequence(A,E)
            if verbosity: print('>random_sequence_%i\n%s\n' % (i,seq))
            f.write('>random_sequence_%i\n%s\n' % (i,seq))

    #####################
    #  END CODING HERE  #
    #####################


if __name__ == "__main__":
    main()
