import argparse
import sys
from pathlib import Path
from.workflow_runner import run_workflow
from .operations.fasta_stats import (
    write_fasta_stats,
    write_assembly_stats
)

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
    
    return 1