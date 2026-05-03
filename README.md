# IBSY Assignment 4 - Sequence Alignment

## Report (PDF)

The final report is available here:

- [Open report](./Latex/Document.pdf)
- [Download report PDF](./Latex/Document.pdf?raw=1)

## Project Overview

This repository contains the implementation and evaluation of sequence alignment methods for IBSY Assignment 4.

Implemented methods:

- Pairwise: Needleman-Wunsch, Smith-Waterman, Gotoh, Banded DP, Hirschberg
- Heuristics: BLAST-style Seed-and-Extend, Minimizer-Greedy
- MSA: Progressive Alignment, Iterative Refinement, Profile-HMM (Viterbi)

## Structure

- [Implementation](Implementation) - Python code and tests
- [Datasets](Datasets) - input data
- [Instructions](Instructions) - assignment description
- [Latex](Latex) - LaTeX sources and final PDF

## Run

Direct execution (script is executable):

```bash
./Implementation/main.py
```

Alternative with Python:

```bash
python Implementation/main.py
```

## Tests

```bash
python -m unittest discover -s Implementation/tests -v
```

## Rebuild PDF

```bash
cd Latex
latexmk -pdf Document.tex
```
