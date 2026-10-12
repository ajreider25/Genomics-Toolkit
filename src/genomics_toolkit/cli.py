import argparse
from.workflow_runner import run_workflow

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
        type=str,
        default=None,
        help="Path to the input FASTA file."
    )

    run_parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Path for the output CSV report."
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
    
    return 1