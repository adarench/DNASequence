Gene Sequencing Alignment Project

Project Overview

This project implements two versions of the Needleman-Wunsch algorithm for gene sequence alignment:

Unrestricted Needleman-Wunsch Algorithm - A classic dynamic programming approach that performs global alignment between two sequences.
Banded Needleman-Wunsch Algorithm - An optimized version that limits computations to a band around the main diagonal, enhancing performance for highly similar sequences.
This project demonstrates practical applications of dynamic programming, performance optimization, and algorithmic problem-solving for bioinformatics.

Key Features

Unrestricted Alignment - Computes the optimal alignment score and sequence alignment using an O(nm) complexity approach.
Banded Alignment - Reduces the computational complexity to O(kn) by restricting the alignment to a narrow band, making it feasible for long sequences with limited deviation from the main diagonal.
Dynamic Programming Implementation - Utilizes a matrix-based approach to store alignment scores and backtracking pointers.
Technologies Used

Python - Core programming language for algorithm implementation.
Git - Version control for tracking changes and managing repository.
Pytest - Testing framework for validating functionality and performance requirements.
Getting Started

Prerequisites
Ensure you have the following installed:

Python 3.x
Git
To install necessary packages, use:

pip install -r requirements.txt
Running the Code
Clone the repository:
git clone https://github.com/your-username/your-repo-name.git
Navigate to the project directory:
cd your-repo-name
Run the main alignment script with example sequences:
python main.py ATCTGG ATGGG
Run the tests:
pytest test_alignment.py
Usage

The primary script, main.py, allows users to input two sequences for alignment. The program outputs the alignment score and the optimal aligned sequences. To change between unrestricted and banded alignment, specify the banded_width parameter in the align function.

Project Structure

alignment.py - Contains the main implementation of the Needleman-Wunsch and banded alignment algorithms.
main.py - Interface to run alignments between two sequences and print results.
test_alignment.py - Test suite to validate correctness and performance.
README.md - Project description and usage information.
Performance Comparison

The unrestricted algorithm offers exact global alignment but has O(nm) time complexity. In contrast, the banded algorithm provides approximate alignment within a limited bandwidth k with an O(kn) complexity, making it more efficient for long, similar sequences.

Algorithm	                       Time Complexity	Space Complexity	Suitable For
Unrestricted Needleman-Wunsch	       O(nm)	        O(nm)	       Short to medium-length sequences
Banded Needleman-Wunsch	               O(kn)	        O(kn)	       Long sequences with minimal deviations

Applications

This project has practical applications in bioinformatics, particularly in:

Gene Sequencing - Comparing genetic material to assess evolutionary relationships or identify gene similarities.
Disease Research - Aligning viral sequences to study mutations and similarities across strains.

Future Improvements

Heuristic Optimization - Implementing additional heuristics for faster approximate alignment.
Graphical User Interface - Adding a GUI to visualize alignment results.
Alternative Scoring Schemes - Experimenting with different scoring penalties to observe their effect on alignment results.