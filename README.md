# Genomics Toolkit

A modular, extensible Python toolkit for genomic data processing, analysis, and reproducible bioinformatics workflows.

## Overview

**Genomics Toolkit** is an open-source project focused on developing computational tools that help researchers process genomic data and investigate biological questions.

The project begins with foundational operations for common bioinformatics file formats and aims to expand into integrated analysis workflows, comparative genomics, and machine-learning applications.

The long-term goal is to provide both **ready-to-use tools for common biological analyses** and **flexible components for building more comprehensive research workflows**.

### Design Principles

- **Biological relevance:** Develop functionality that supports meaningful biological analysis.
- **Modularity:** Keep individual operations reusable and independently testable.
- **Reproducibility:** Support consistent, documented, and repeatable analyses.
- **Extensibility:** Allow new file formats, operations, and analytical methods to be integrated.
- **Scalability:** Design with larger datasets, parallel execution, and comprehensive workflows in mind.
- **Accessibility:** Provide a straightforward command-line interface alongside reusable Python functions.

## Features

### FASTA Operations

The toolkit provides foundational functionality for working with FASTA sequence files.

| Command | Description |
|---|---|
| `fasta-stats` | Calculate sequence statistics and assembly-level metrics. |
| `fasta-validate` | Validate FASTA structure and nucleotide sequences. |
| `fasta-filter` | Filter sequences by length and identifiers. |
| `fasta-revcomp` | Generate reverse-complement sequences. |
| `fasta-extract` | Extract sequences using a file containing sequence identifiers. |
| `fasta-split` | Split FASTA datasets into smaller files by record count. |
| `fasta-merge` | Combine multiple FASTA files into one dataset. |

FASTA operations are designed to support both uncompressed and gzip-compressed data where applicable.

### FASTQ Operations - in development

FASTQ functionality is under active development.

| Command | Description | Status |
|---|---|---|
| `fastq-validate` | Validate FASTQ records, nucleotide sequences, and quality-score structure. | Implemented |
| `fastq-stats` | Calculate sequencing-read and quality statistics. | Planned |
| `fastq-filter` | Filter reads by length, quality, and other criteria. | Planned |
| `fastq-trim` | Trim low-quality bases from sequencing reads. | Planned |
| `fastq-sample` | Sample reads from FASTQ datasets. | Planned |
| `fastq-convert` | Convert FASTQ sequences to FASTA. | Planned |
| `fastq-merge` | Combine FASTQ datasets. | Planned |

## Installation

### Prerequisites

- Python 3.14 or a compatible version supported by the project
- Git
- [uv](https://docs.astral.sh/uv/), the Python package and project manager

Clone the repository:

```bash
git clone https://github.com/<your-username>/Genomics-Toolkit.git
cd Genomics-Toolkit
```

Install project dependencies:

```bash
uv sync
```

Verify that the CLI is available:

```bash
uv run genomics-toolkit --help
```

Optional workflow dependencies can be installed with:

```bash
uv sync --group workflow
```

## Usage

The toolkit provides a command-line interface for executing individual genomic operations.

### FASTA Statistics

```bash
uv run genomics-toolkit fasta-stats \
    --input data/sample.fasta
```

### FASTA Validation

```bash
uv run genomics-toolkit fasta-validate \
    --input data/sample.fasta \
    --strict \
    --report results/fasta_validation.json
```

### FASTA Filtering

```bash
uv run genomics-toolkit fasta-filter \
    --input data/sample.fasta \
    --output results/filtered.fasta \
    --min-length 100
```

### FASTA Reverse Complement

```bash
uv run genomics-toolkit fasta-revcomp \
    --input data/sample.fasta \
    --output results/revcomp.fasta
```

### FASTA Extraction

```bash
uv run genomics-toolkit fasta-extract \
    --input data/sample.fasta \
    --ids-file data/selected_ids.txt \
    --output results/extracted.fasta
```

### FASTA Splitting

```bash
uv run genomics-toolkit fasta-split \
    --input data/sample.fasta \
    --output-dir results/fasta_chunks \
    --records-per-file 1000
```

### FASTA Merging

```bash
uv run genomics-toolkit fasta-merge \
    --inputs results/fasta_chunks/part_001.fasta \
             results/fasta_chunks/part_002.fasta \
    --output results/merged.fasta
```

### FASTQ Validation

```bash
uv run genomics-toolkit fastq-validate \
    --input data/sample.fastq \
    --report results/fastq_validation.json
```

For command-specific options, use:

```bash
uv run genomics-toolkit <command> --help
```

## Project Structure

```text
Genomics-Toolkit/
├── src/
│   └── genomics_toolkit/
│       ├── cli.py
│       ├── workflow_runner.py
│       └── operations/
│           ├── fasta/
│           │   ├── fasta_stats.py
│           │   ├── fasta_validate.py
│           │   ├── fasta_filter.py
│           │   ├── fasta_revcomp.py
│           │   ├── fasta_extract.py
│           │   ├── fasta_split.py
│           │   └── fasta_merge.py
│           └── fastq/
│               └── fastq_validate.py
├── tests/
│   ├── fasta/
│   └── fastq/
├── workflow/
├── data/
├── pyproject.toml
└── README.md
```

The project separates reusable operations from command-line handling and workflow orchestration.

## Testing

Automated testing is performed using `pytest`.

Run the complete test suite:

```bash
uv run pytest -v
```

Run FASTA-specific tests:

```bash
uv run pytest tests/fasta/ -v
```

Run FASTQ-specific tests:

```bash
uv run pytest tests/fastq/ -v
```

New operations should include tests covering expected behavior, malformed input, edge cases, and file-handling safety.

## Development Workflow

Development follows a feature-branch workflow.

Individual implementations are developed and tested on separate branches.

For larger development phases, such as FASTQ, individual feature branches are merged into a format-specific integration branch.

For example:

```text
main
└── feature/fastq
    ├── feature/fastq-validate
    ├── feature/fastq-stats
    ├── feature/fastq-filter
    └── ...
```

Once a format's foundational functionality is complete and tested, its integration branch can be merged into `main`.

This approach helps keep the main branch stable while allowing incremental development.

## Development Roadmap

### Phase 1: Foundational Genomics Operations

Establish reliable processing and analysis functionality for widely used genomics formats.

| Format | Focus |
|---|---|
| FASTA | Sequence processing and assembly statistics |
| FASTQ | Sequencing-read processing and quality control |
| VCF/BCF | Genetic variation and genotype analysis |
| SAM/BAM/CRAM | Sequencing alignments and coverage |
| BED | Genomic interval operations |
| GFF3/GTF | Gene and genomic feature annotations |

### Phase 2: Biological Analysis Workflows

Develop configurable, reproducible pipelines for common biological applications, including:

- Genome assembly quality assessment
- Sequencing-read quality control
- Alignment quality assessment
- Variant analysis
- Comparative genomics
- Genomic annotation analysis

The toolkit includes initial workflow infrastructure using Snakemake. More comprehensive workflows are planned for future development.

### Phase 3: Advanced Analysis and Machine Learning

Extend the toolkit with capabilities for:

- Genomic feature extraction and engineering
- Sequence classification
- Variant effect prediction
- Comparative genomic analysis
- Biological phenotype prediction
- Reproducible machine-learning experiments

Machine learning will complement established bioinformatics methods rather than replace them.

## Project Status

**Active development — foundational genomics operations.**

The project is currently focused on completing its core file-processing and analysis functionality before expanding into comprehensive biological workflows and machine-learning applications.

Interfaces, supported formats, and functionality may evolve as development progresses.

## Contributing

Contributions, bug reports, feature suggestions, and improvements are welcome.

When contributing:

1. Create a branch for your proposed change.
2. Keep implementations modular and focused.
3. Include or update relevant tests.
4. Run the test suite before submitting changes.
5. Open a pull request describing the change and its purpose.

## License

License information will be added once a license has been selected.