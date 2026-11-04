# Genomics Toolkit

**An AI-guided, hands-on bioinformatics learning project for understanding genomic data formats, developing computational tools, and exploring how bioinformatics is used to solve biological problems.**

## Overview

Genomics Toolkit is a Python-based project focused on learning bioinformatics through practical software development.

The project explores the major file formats used in genomics, the biological information they represent, and the computational operations used to analyze them.

Rather than simply using existing bioinformatics software, I am building my own implementations of common operations to develop a deeper understanding of how genomic data is structured, processed, validated, and interpreted.

Development is guided by AI-assisted learning, with an emphasis on understanding the underlying concepts, reviewing implementations, writing tests, and gradually building more sophisticated tools.

The long-term goal is to connect these foundational operations into reproducible workflows capable of addressing real biological questions.

## Why I'm Building This

This project began as an effort to create a collection of genomics utilities, but it has evolved into a personal, AI-guided bioinformatics education project.

My goal is to learn the different data formats used throughout genomics, understand what each format represents, and explore how bioinformaticians use these data to investigate biological problems.

I want to go beyond memorizing file specifications or learning commands. By implementing the operations myself, I can better understand what happens to genomic data at each stage of an analysis, why particular processing steps are necessary, and how individual tools fit into larger biological workflows.

AI serves as a learning and development partner throughout this process, helping me explore concepts, design implementations, troubleshoot problems, and evaluate different approaches. The objective is not simply to generate working code, but to develop my own understanding of the biology, data structures, algorithms, and software engineering practices involved.

This repository documents that learning process as the toolkit grows.

## Learning Objectives

Through this project, I aim to:

- Understand the major file formats used in genomics and the biological information they contain.
- Learn how common bioinformatics operations work by implementing them in Python.
- Explore how genomic data moves through different stages of an analysis.
- Develop skills in Python programming, command-line interfaces, automated testing, and version control.
- Understand the biological reasoning behind data processing and analysis decisions.
- Build reusable and reliable tools that can eventually be combined into complete bioinformatics workflows.
- Explore more advanced approaches, including workflow automation, scalable analysis, and machine learning, as my understanding develops.

## Supported Data Formats

Development is organized around the major data formats encountered in genomics.

| Format | Biological Information | Development Status |
|---|---|---|
| FASTA | DNA, RNA, and protein sequences | Foundational operations implemented |
| FASTQ | Sequencing reads and base quality scores | Foundational operations implemented |
| VCF/BCF | Genetic variants and genotype information | Planned next |
| SAM/BAM/CRAM | Sequence alignments against reference genomes | Planned |
| BED | Genomic intervals and regions | Planned |
| GFF3/GTF | Gene structures and genomic annotations | Planned |

Each format is approached as both a software development milestone and a learning opportunity.

## Implemented Features

### FASTA — Sequence Data

FASTA is commonly used to represent nucleotide and protein sequences.

Implemented operations:

| Command | Purpose |
|---|---|
| `fasta-stats` | Calculate sequence and assembly statistics |
| `fasta-validate` | Validate FASTA structure and sequence content |
| `fasta-filter` | Select sequences using defined criteria |
| `fasta-revcomp` | Generate reverse-complement sequences |
| `fasta-extract` | Extract sequences by identifier |
| `fasta-split` | Divide FASTA files into smaller files |
| `fasta-merge` | Combine multiple FASTA files |

**Learning focus:** Sequence representation, parsing, sequence composition, genome assembly statistics, and basic sequence manipulation.

### FASTQ — Sequencing Read Data

FASTQ represents sequencing reads together with their corresponding base quality scores.

Implemented operations:

| Command | Purpose |
|---|---|
| `fastq-validate` | Validate FASTQ records and quality scores |
| `fastq-stats` | Calculate sequencing read and quality statistics |
| `fastq-filter` | Filter reads using length, quality, and ambiguous-base criteria |
| `fastq-trim` | Trim low-quality bases from read ends |
| `fastq-sample` | Randomly sample sequencing reads |
| `fastq-convert` | Convert FASTQ records to FASTA |
| `fastq-merge` | Combine multiple FASTQ files |
| `fastq-pair-validate` | Check paired-end FASTQ synchronization |

**Learning focus:** Sequencing quality scores, read preprocessing, quality control, paired-end sequencing, and preparation of data for downstream analysis.

Additional paired-end filtering and sampling operations are potential future enhancements.

## Installation

The project uses Python and [uv](https://docs.astral.sh/uv/) for dependency management.

Clone the repository:

```bash
git clone <repository-url>
cd Genomics-Toolkit
```

Install dependencies:

```bash
uv sync
```

Display the command-line interface:

```bash
uv run genomics-toolkit --help
```

## Example Usage

### Validate a FASTA file

```bash
uv run genomics-toolkit fasta-validate \
  --input data/sample.fasta
```

### Calculate FASTA statistics

```bash
uv run genomics-toolkit fasta-stats \
  --input data/sample.fasta \
  --output results/fasta_stats.json
```

### Validate FASTQ sequencing reads

```bash
uv run genomics-toolkit fastq-validate \
  --input data/sample.fastq
```

### Trim low-quality reads

```bash
uv run genomics-toolkit fastq-trim \
  --input data/sample.fastq \
  --output results/trimmed.fastq \
  --min-quality 20 \
  --min-length 4
```

### Filter sequencing reads

```bash
uv run genomics-toolkit fastq-filter \
  --input results/trimmed.fastq \
  --output results/filtered.fastq \
  --min-length 4 \
  --max-n-percent 20
```

### Calculate FASTQ statistics

```bash
uv run genomics-toolkit fastq-stats \
  --input results/filtered.fastq \
  --output results/filtered_stats.json
```

### Validate paired-end sequencing reads

```bash
uv run genomics-toolkit fastq-pair-validate \
  --r1 data/sample_R1.fastq \
  --r2 data/sample_R2.fastq
```

These examples demonstrate how individual operations can be combined into simple data-processing workflows.

## Development Approach

The project follows an incremental, format-by-format development strategy.

1. **Understand the format:** Study its structure, biological meaning, and common applications.
2. **Implement foundational operations:** Build reusable Python functions for validation, statistics, filtering, and manipulation.
3. **Expose CLI commands:** Make operations accessible through a consistent command-line interface.
4. **Write automated tests:** Verify expected behavior, edge cases, and error handling using pytest.
5. **Integrate operations:** Explore how individual tools work together to solve larger analysis problems.
6. **Expand into workflows:** Gradually introduce reproducible, automated pipelines that connect multiple formats and analysis stages.

This approach prioritizes understanding and reliability before increasing complexity.

## Testing

Automated tests are written using pytest.

Run the complete test suite:

```bash
uv run pytest -v
```

Run tests for a particular format:

```bash
uv run pytest tests/fasta/ -v
uv run pytest tests/fastq/ -v
```

Tests help ensure that operations behave consistently as new functionality is introduced.

## Roadmap

### Phase 1 — Foundational Genomic Data Formats

- [x] FASTA sequence operations
- [x] FASTQ sequencing read operations
- [ ] VCF/BCF variant operations
- [ ] SAM/BAM/CRAM alignment operations
- [ ] BED genomic interval operations
- [ ] GFF3/GTF annotation operations

### Phase 2 — Integration and Biological Workflows

- Develop reusable operations that work across genomic data formats.
- Create standard analysis workflows for common bioinformatics tasks.
- Integrate workflow management and reproducible execution.
- Explore realistic biological questions using representative datasets.

### Phase 3 — Advanced Capabilities

- Expand workflow functionality and scalability.
- Explore additional genomic data formats.
- Investigate computational methods for more advanced analyses.
- Explore machine learning applications where biologically appropriate.

The roadmap is expected to evolve alongside the learning process.

## Project Philosophy

The guiding principle of this project is **learning bioinformatics by building bioinformatics tools**.

Working software is important, but understanding the reasoning behind that software is equally important.

Each feature should contribute to a broader understanding of how genomic data is generated, represented, transformed, and analyzed.

Ultimately, the goal is not simply to accumulate commands. It is to develop the knowledge and practical skills needed to connect computational methods with meaningful biological questions.

## AI-Assisted Development

This project is developed with assistance from AI tools, including ChatGPT.

AI is used to support conceptual learning, software design, code development, debugging, testing, and technical discussion.

The project is intended to be an active learning process rather than a collection of automatically generated solutions. Understanding, testing, and critically evaluating implementations remain central to its purpose.

## Project Status

This is an evolving educational project. Implementations may be simplified compared with mature production bioinformatics software and should not be assumed to be suitable for clinical, diagnostic, or other high-stakes applications without independent validation.