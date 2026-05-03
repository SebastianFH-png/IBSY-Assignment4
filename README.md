# IBSY Assignment 4 - Sequence Alignment

## Report (PDF)

The final report is available here:

- [Open report](./Latex/Document.pdf)
- [Download report PDF](./Latex/Document.pdf?raw=1)

## Project Overview

This repository contains the implementation and evaluation of sequence alignment methods for IBSY Assignment 4.

Implemented methods:

- Needleman-Wunsch (Global Alignment)
- Smith-Waterman (Local Alignment)
- Gotoh (Global Alignment with Affine Gaps)
- Basic Seed-and-Extend (Heuristic / BLAST-style)

## Structure

- [Implementation](Implementation) - Python code
- [Datasets](Datasets) - input data (fasta format)
- [Instructions](Instructions) - assignment description
- [Latex](Latex) - LaTeX sources for the report

## Run

To run the experiments:

```bash
python Implementation/Main.py
```

## Rebuild PDF

```bash
cd Latex
latexmk -pdf Document.tex
```
