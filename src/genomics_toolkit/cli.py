import argparse
import sys
from pathlib import Path
from.workflow_runner import run_workflow
from .operations.fasta_stats import (
    write_fasta_stats,
    write_assembly_stats
)
from .operations.fasta_validate import (
    validate_fasta,
    write_validation_report
)
from .operations.fasta_filter import filter_fasta
from .operations.fasta_revcomp import write_reverse_complement


def build_parser() -> argparse.ArgumentParser:
    """Create the command-line argument parser."""

    parser = argparse.ArgumentParser(
        prog="genomics",
        description="A modular toolkit for genomic analysis."
    )

    subparsers = parser.add_subparsers(
        dest="command",
        required=True
    )

    ## Pipeline execution command
    run_parser = subparsers.add_parser(
        "run",
        help="Run a genomics workflow."
    )

    run_parser.add_argument(
        "--cores",
        type=int,
        default=1,
        help="Number of CPU cores to use (def: 1)"
    )

    run_parser.add_argument(
        "--snakefile",
        type=str,
        default="workflow/Snakefile",
        help="Path to the Snakefile."
    )

    run_parser.add_argument(
        "--input",
        default=None,
        help="Input FASTA file.",
    )

    run_parser.add_argument(
        "--output",
        default=None,
        help="Output CSV file.",
    )

    fasta_parser = subparsers.add_parser(
        "fasta-stats",
        help="Calculate statistics for a FASTA file."
    )

    fasta_parser.add_argument(
        "--input",
        required=True,
        help="Path to the input FASTA file."
    )

    fasta_parser.add_argument(
        "--output",
        required=True,
        help="Path for the output CSV report."
    )

    fasta_parser.add_argument(
        "--summary",
        default=None,
        help="Optional path for assembly summary CSV."
    )

    validate_parser = subparsers.add_parser(
        "fasta-validate",
        help="Validate DNA FASTA files."
    )

    validate_parser.add_argument(
        "--input",
        required=True,
        help="Path to the FASTA file."
    )

    validate_parser.add_argument(
        "--strict",
        action="store_true",
        help="Treat empty sequences and empty files as errors."
    )

    validate_parser.add_argument(
        "--report",
        default=None,
        help="Optional path for a JSON validation report."
    )

    filter_parser = subparsers.add_parser(
        "fasta-filter",
        help="Filter FASTA records by length or sequence id."
    )

    filter_parser.add_argument(
        "--input",
        required=True,
        help="Input FASTA file."
    )

    filter_parser.add_argument(
        "--output",
        required=True,
        help="Output FASTA file"
    )

    filter_parser.add_argument(
        "--min-length",
        type=int,
        default=None,
        help="Minimum sequence length (inclusive)."
    )

    filter_parser.add_argument(
        "--max-length",
        type=int,
        default=None,
        help="Maximum sequence length (inclusive)."
    )

    filter_parser.add_argument(
        "--ids",
        default=None,
        help="Comma-seperated sequence IDs to retain."
    )

    revcomp_parser = subparsers.add_parser(
        "fasta-revcomp",
        help="Generate reverse-complement DNA FASTA sequences."
    )

    revcomp_parser.add_argument(
        "--input",
        required=True,
        help="Input FASTA or FASTA.gz file"
    )

    revcomp_parser.add_argument(
        "--output",
        required=True,
        help="Output FASTA or FASTA.gz file"
    )

    return parser

def main(argv: list[str] | None = None) -> int:
    """Execute the requested CLI command."""

    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "run":
        if args.cores < 1:
            parser.error("--cores must be at least 1")
    
        return run_workflow(
            cores=args.cores,
            snakefile=args.snakefile,
            input_path=args.input,
            output_path=args.output
        )

    if args.command == "fasta-stats":
        try:
            if args.summary:
                if Path(args.output).resolve() == Path(args.summary).resolve():
                    print(
                        "Error: --output and --summary must be different files.",
                        file=sys.stderr,
                    )
                    return 1
                
            write_fasta_stats(
                input_path=args.input,
                output_path=args.output
            )

            if args.summary:
                write_assembly_stats(
                    input_path=args.input,
                    output_path=args.summary
                )

        except (OSError, ValueError) as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1
    
        print(f"FASTA statistics written to: {args.output}")

        if args.summary:
            print(f"Assembly summary written to: {args.summary}")

        return 0

    if args.command == "fasta-validate":
        try:
            if args.report:
                if Path(args.input).resolve() == Path(args.report).resolve():
                    print(
                        "Error: Input FASTA and JSON report must be different files.",
                        file=sys.stderr
                    )
                    return 1
                
            report = validate_fasta(
                input_path=args.input,
                strict=args.strict
            )

            if args.report:
                write_validation_report(
                    report=report,
                    output_path=args.report,
                    input_path=args.input)
                
        except (OSError, ValueError, EOFError) as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1

        print(f"Sequences checked: {report['sequence_count']}")
        print(f"Errors: {report['error_count']}")
        print(f"Warnings: {report['warning_count']}")

        for error in report["errors"]:
            print(f"ERROR: {error}")

        for warning in report["warnings"]:
            print(f"WARNING: {warning}")

        if report["valid"]:
            print("FASTA validation passed.")
            return 0

        print("FASTA validation failed.")
        return 1

    if args.command == "fasta-filter":
        try:
            selected_ids = None

            if args.ids is not None:
                selected_ids = {
                    item.strip()
                    for item in args.ids.split(",")
                    if item.strip()
                }

                if not selected_ids:
                    raise ValueError(
                        "--ids must contain at least one sequence ID."
                    )

            result = filter_fasta(
                input_path=args.input,
                output_path=args.output,
                min_length=args.min_length,
                max_length=args.max_length,
                ids=selected_ids
            )

        except (OSError, ValueError, EOFError) as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1

        print(f"Sequences processed: {result['total_records']}")
        print(f"Sequences retained: {result['retained_records']}")
        print(f"Sequences removed: {result['removed_records']}")

        return 0

    if args.command == "fasta-revcomp":
        try:
            result = write_reverse_complement(
                input_path=args.input,
                output_path=args.output
            )

        except (OSError, ValueError, EOFError) as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1

        print(f"Sequences processed: {result['sequence_count']}")
        print(f"Total bases: {result['total_bases']}")

        return 0
    
    return 1