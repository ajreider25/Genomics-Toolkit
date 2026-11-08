import argparse
import sys
from pathlib import Path
from.workflow_runner import run_workflow
from .operations.fasta.fasta_stats import (
    write_fasta_stats,
    write_assembly_stats
)
from .operations.fasta.fasta_validate import (
    validate_fasta,
    write_validation_report
)
from .operations.fasta.fasta_filter import filter_fasta
from .operations.fasta.fasta_revcomp import write_reverse_complement
from .operations.fasta.fasta_extract import extract_fasta
from .operations.fasta.fasta_split import split_fasta
from .operations.fasta.fasta_merge import merge_fasta

from .operations.fastq.fastq_validate import (
    validate_fastq,
    write_validation_report as write_fastq_validation_report
)
from .operations.fastq.fastq_stats import (
    calculate_fastq_stats,
    write_fastq_stats
)
from .operations.fastq.fastq_filter import filter_fastq
from .operations.fastq.fastq_trim import trim_fastq
from .operations.fastq.fastq_sample import sample_fastq
from .operations.fastq.fastq_convert import fastq_to_fasta
from .operations.fastq.fastq_merge import merge_fastq
from .operations.fastq.fastq_pair_validate import validate_fastq_pairs

from .operations.vcf.vcf_validate import validate_vcf
from .operations.vcf.vcf_stats import (
    calculate_vcf_stats,
    write_vcf_stats
)
from .operations.vcf.vcf_filter import filter_vcf
from .operations.vcf.vcf_extract import extract_vcf

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



    stats_parser = subparsers.add_parser(
        "fasta-stats",
        help="Calculate statistics for a FASTA file."
    )

    stats_parser.add_argument(
        "--input",
        required=True,
        help="Path to the input FASTA file."
    )

    stats_parser.add_argument(
        "--output",
        required=True,
        help="Path for the output CSV report."
    )

    stats_parser.add_argument(
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



    extract_parser = subparsers.add_parser(
        "fasta-extract",
        help="Extract FASTA sequences using an ID selection file."
    )

    extract_parser.add_argument(
        "--input",
        required=True,
        help="Text FASTA or FASTA.gz file."
    )

    extract_parser.add_argument(
        "--ids-file",
        required=True,
        help="Text file containong one sequence ID per line."
    )

    extract_parser.add_argument(
        "--output",
        required=True,
        help="Output FASTA or FASTA.gz file."
    )



    split_parser = subparsers.add_parser(
        "fasta-split",
        help="Split a FASTA into smaller files."
    )

    split_parser.add_argument(
        "--input",
        required=True,
        help="Input FASTA or FASTA.gz file."
    )

    split_parser.add_argument(
        "--output-dir",
        required=True,
        help="New directory for output FASTA chunks."
    )

    split_parser.add_argument(
        "--records-per-file",
        type=int,
        default=1000,
        help="Maximum number of records per file (default: 1000)."
    )

    split_parser.add_argument(
        "--gzip",
        action="store_true",
        help="Compress output files using gzip."
    )


    merge_parser = subparsers.add_parser(
        "fasta-merge",
        help="Merge multiple FASTA files.",
    )

    merge_parser.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="Input FASTA files in merge order.",
    )

    merge_parser.add_argument(
        "--output",
        required=True,
        help="Output FASTA or FASTA.gz file.",
    )

    merge_parser.add_argument(
        "--allow-duplicates",
        action="store_true",
        help="Allow duplicate sequence IDs.",
    )


    fastq_validate_parser = subparsers.add_parser(
        "fastq-validate",
        help="Validate FASTQ files.",
    )

    fastq_validate_parser.add_argument(
        "--input",
        required=True,
        help="Input FASTQ or FASTQ.gz file.",
    )

    fastq_validate_parser.add_argument(
        "--strict",
        action="store_true",
        help="Treat empty sequences and empty files as errors.",
    )

    fastq_validate_parser.add_argument(
        "--report",
        help="Optional JSON validation report path.",
    )


    
    fastq_stats_parser = subparsers.add_parser(
        "fastq-stats",
        help="Calculate FASTQ sequencing quality statistics.",
    )

    fastq_stats_parser.add_argument(
        "--input",
        required=True,
        help="Input FASTQ or FASTQ.gz file.",
    )

    fastq_stats_parser.add_argument(
        "--output",
        help="Optional JSON statistics output file.",
    )


    
    fastq_filter_parser = subparsers.add_parser(
        "fastq-filter",
        help="Filter FASTQ reads by length, quality, and N content.",
    )

    fastq_filter_parser.add_argument("--input", required=True)
    fastq_filter_parser.add_argument("--output", required=True)
    fastq_filter_parser.add_argument("--min-length", type=int)
    fastq_filter_parser.add_argument("--max-length", type=int)
    fastq_filter_parser.add_argument("--min-quality", type=float)
    fastq_filter_parser.add_argument("--max-n-percent", type=float)



    fastq_trim_parser = subparsers.add_parser(
        "fastq-trim",
        help="Trim low-quality bases from the 3' ends of FASTQ reads.",
    )

    fastq_trim_parser.add_argument(
        "--input",
        required=True,
        help="Input FASTQ or FASTQ.gz file.",
    )

    fastq_trim_parser.add_argument(
        "--output",
        required=True,
        help="Output FASTQ or FASTQ.gz file.",
    )

    fastq_trim_parser.add_argument(
        "--min-quality",
        type=int,
        default=20,
        help="Minimum Phred quality to retain at the 3' end (default: 20).",
    )

    fastq_trim_parser.add_argument(
        "--min-length",
        type=int,
        default=0,
        help="Discard reads shorter than this after trimming (default: 0).",
    )


    fastq_sample_parser = subparsers.add_parser(
        "fastq-sample",
        help="Randomly sample FASTQ reads.",
    )

    fastq_sample_parser.add_argument(
        "--input",
        required=True,
        help="Input FASTQ or FASTQ.gz file.",
    )

    fastq_sample_parser.add_argument(
        "--output",
        required=True,
        help="Output FASTQ or FASTQ.gz file.",
    )

    sample_group = fastq_sample_parser.add_mutually_exclusive_group(
        required=True
    )

    sample_group.add_argument(
        "--count",
        type=int,
        help="Exact number of reads to sample (or all reads if fewer exist).",
    )

    sample_group.add_argument(
        "--fraction",
        type=float,
        help="Probability of selecting each read, from 0 to 1.",
    )

    fastq_sample_parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducible sampling (default: 42).",
    )


    fastq_convert_parser = subparsers.add_parser(
        "fastq-convert",
        help="Convert FASTQ reads to FASTA sequences.",
    )

    fastq_convert_parser.add_argument(
        "--input",
        required=True,
        help="Input FASTQ or FASTQ.gz file.",
    )

    fastq_convert_parser.add_argument(
        "--output",
        required=True,
        help="Output FASTA or FASTA.gz file.",
    )

    fastq_convert_parser.add_argument(
        "--line-width",
        type=int,
        default=80,
        help="FASTA sequence line width (default: 80).",
    )

 
    fastq_merge_parser = subparsers.add_parser(
        "fastq-merge",
        help="Merge multiple FASTQ files.",
    )

    fastq_merge_parser.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="Input FASTQ or FASTQ.gz files in merge order.",
    )

    fastq_merge_parser.add_argument(
        "--output",
        required=True,
        help="Output FASTQ or FASTQ.gz file.",
    )


    fastq_pair_validate_parser = subparsers.add_parser(
        "fastq-pair-validate",
        help="Validate synchronized paired-end FASTQ files.",
    )

    fastq_pair_validate_parser.add_argument(
        "--r1",
        required=True,
        help="R1 FASTQ or FASTQ.gz file.",
    )

    fastq_pair_validate_parser.add_argument(
        "--r2",
        required=True,
        help="R2 FASTQ or FASTQ.gz file.",
    )


    vcf_validate_parser = subparsers.add_parser(
        "vcf-validate",
        help="Validate the structure of a VCF file.",
    )

    vcf_validate_parser.add_argument(
        "--input",
        required=True,
        help="Input VCF file to validate.",
    )


    vcf_stats_parser = subparsers.add_parser(
        "vcf-stats",
        help="Calculate summary statistics for a VCF file.",
    )

    vcf_stats_parser.add_argument(
        "--input",
        required=True,
        help="Input VCF file.",
    )

    vcf_stats_parser.add_argument(
        "--output",
        required=True,
        help="Output JSON statistics file.",
    )



    vcf_filter_parser = subparsers.add_parser(
        "vcf-filter",
        help="Filter VCF records by quality, depth, and variant type.",
    )

    vcf_filter_parser.add_argument(
        "--input",
        required=True,
        help="Input VCF file.",
    )

    vcf_filter_parser.add_argument(
        "--output",
        required=True,
        help="Output filtered VCF file.",
    )

    vcf_filter_parser.add_argument(
        "--min-qual",
        type=float,
        help="Minimum variant QUAL score.",
    )

    vcf_filter_parser.add_argument(
        "--pass-only",
        action="store_true",
        help="Keep only records with FILTER=PASS.",
    )

    vcf_filter_parser.add_argument(
        "--chrom",
        help="Keep variants on this chromosome.",
    )

    vcf_filter_parser.add_argument(
        "--min-depth",
        type=int,
        help="Minimum site-level INFO/DP.",
    )

    vcf_filter_parser.add_argument(
        "--variant-type",
        choices=[
            "snv",
            "mnv",
            "insertion",
            "deletion",
            "complex",
        ],
        help="Keep records containing this variant type.",
    )

 
    vcf_extract_parser = subparsers.add_parser(
        "vcf-extract",
        help="Extract VCF records by region or variant ID.",
    )

    vcf_extract_parser.add_argument(
        "--input",
        required=True,
        help="Input VCF file.",
    )

    vcf_extract_parser.add_argument(
        "--output",
        required=True,
        help="Output VCF file.",
    )

    vcf_extract_parser.add_argument(
        "--chrom",
        help="Chromosome to extract.",
    )

    vcf_extract_parser.add_argument(
        "--start",
        type=int,
        help="Start position (1-based, inclusive).",
    )

    vcf_extract_parser.add_argument(
        "--end",
        type=int,
        help="End position (1-based, inclusive).",
    )

    vcf_extract_parser.add_argument(
        "--id",
        dest="variant_id",
        help="Variant identifier to extract.",
    )

    vcf_extract_parser.add_argument(
        "--ids-file",
        help="Text file containing variant IDs.",
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

    if args.command == "fasta-extract":
        try:
            result = extract_fasta(
                input_path=args.input,
                ids_file=args.ids_file,
                output_path=args.output
            )

        except (OSError, ValueError, EOFError) as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1

        print(f"Sequences processed: {result['total_records']}")
        print(f"Sequences extracted: {result['extracted_records']}")
        print(f"Unique IDs requested: {result['requested_ids']}")
        print(f"Unique IDs found: {result['found_ids']}")

        if result["missing_ids"]:
            print(
                "WARNING: IDs not found: "
                + ", ".join(result["missing_ids"]),
                file=sys.stderr
            )

        return 0

    if args.command == "fasta-split":
        try:
            result = split_fasta(
                input_path=args.input,
                output_dir=args.output_dir,
                records_per_file=args.records_per_file,
                gzip_output=args.gzip
            )
        except (OSError, ValueError, EOFError) as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1

        print(f"Sequences processed: {result['total_records']}")
        print(f"Output files created: {result['file_count']}")
        print(f"Output directory: {result['output_dir']}")

        return 0

    if args.command == "fasta-merge":
        try:
            result = merge_fasta(
                input_paths=args.inputs,
                output_path=args.output,
            )
        except (OSError, ValueError, UnicodeError) as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1

        print(f"Input files merged: {result['input_files']}")
        print(f"Sequences merged: {result['total_records']}")
        print(f"Total bases: {result['total_bases']}")
        print(f"Duplicate IDs: {result['duplicate_ids']}")
        print(f"Output: {result['output_path']}")

        return 0


    
    if args.command == "fastq-validate":
        try:
            report = validate_fastq(
                input_path=args.input,
                strict=args.strict,
            )

            if args.report:
                write_fastq_validation_report(
                    report,
                    args.report,
                    input_path=args.input,
                )

        except (OSError, ValueError) as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1

        print(f"Valid FASTQ: {report['valid']}")
        print(f"Records: {report['record_count']}")
        print(f"Total bases: {report['total_bases']}")
        print(f"Errors: {report['error_count']}")
        print(f"Warnings: {report['warning_count']}")

        for error in report["errors"]:
            print(f"ERROR: {error}", file=sys.stderr)

        for warning in report["warnings"]:
            print(f"WARNING: {warning}", file=sys.stderr)

        return 0 if report["valid"] else 1

    
    if args.command == "fastq-stats":
        try:
            stats = calculate_fastq_stats(args.input)

            if args.output:
                write_fastq_stats(
                    stats,
                    args.output,
                    input_path=args.input,
                )

        except (OSError, ValueError, UnicodeError) as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1

        print(f"Reads: {stats['read_count']}")
        print(f"Total bases: {stats['total_bases']}")
        print(f"Min read length: {stats['min_read_length']}")
        print(f"Max read length: {stats['max_read_length']}")
        print(f"Mean read length: {stats['mean_read_length']:.2f}")
        print(f"GC content: {stats['gc_percent']:.2f}%")
        print(f"Mean Phred quality: {stats['mean_quality']:.2f}")
        print(f"Q20 bases: {stats['q20_percent']:.2f}%")
        print(f"Q30 bases: {stats['q30_percent']:.2f}%")

        return 0

    
    if args.command == "fastq-filter":
        try:
            result = filter_fastq(
                input_path=args.input,
                output_path=args.output,
                min_length=args.min_length,
                max_length=args.max_length,
                min_quality=args.min_quality,
                max_n_percent=args.max_n_percent,
            )
        except (OSError, ValueError, UnicodeError) as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1

        print(f"Total reads: {result['total_reads']}")
        print(f"Kept reads: {result['kept_reads']}")
        print(f"Filtered reads: {result['filtered_reads']}")
        print(f"Output: {result['output_path']}")

        return 0


    if args.command == "fastq-trim":
        try:
            result = trim_fastq(
                input_path=args.input,
                output_path=args.output,
                min_quality=args.min_quality,
                min_length=args.min_length,
            )

        except (OSError, ValueError, UnicodeError) as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1

        print(f"Total reads: {result['total_reads']}")
        print(f"Kept reads: {result['kept_reads']}")
        print(f"Trimmed reads: {result['trimmed_reads']}")
        print(f"Discarded reads: {result['discarded_reads']}")
        print(f"Bases removed: {result['bases_removed']}")
        print(f"Output: {result['output_path']}")

        return 0

    
    if args.command == "fastq-sample":
        try:
            result = sample_fastq(
                input_path=args.input,
                output_path=args.output,
                count=args.count,
                fraction=args.fraction,
                seed=args.seed,
            )

        except (OSError, ValueError, UnicodeError) as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1

        print(f"Total reads: {result['total_reads']}")
        print(f"Sampled reads: {result['sampled_reads']}")
        print(f"Discarded reads: {result['discarded_reads']}")
        print(f"Method: {result['method']}")
        print(f"Seed: {result['seed']}")
        print(f"Output: {result['output_path']}")

        return 0


    
    if args.command == "fastq-convert":
        try:
            result = fastq_to_fasta(
                input_path=args.input,
                output_path=args.output,
                line_width=args.line_width,
            )

        except (OSError, ValueError, UnicodeError) as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1

        print(f"Converted reads: {result['record_count']}")
        print(f"Total bases: {result['total_bases']}")
        print(f"Output: {result['output_path']}")

        return 0


    if args.command == "fastq-merge":
        try:
            result = merge_fastq(
                input_paths=args.inputs,
                output_path=args.output,
            )
        except (OSError, ValueError, UnicodeError) as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1

        print(f"Input files: {result['input_files']}")
        print(f"Total reads: {result['total_reads']}")
        print(f"Total bases: {result['total_bases']}")
        print(f"Output: {result['output_path']}")
        return 0

    
    if args.command == "fastq-pair-validate":
        result = validate_fastq_pairs(
            r1_path=args.r1,
            r2_path=args.r2,
        )

        print(f"Valid pairs: {result['pair_count']}")
        print(f"R1 bases: {result['r1_bases']}")
        print(f"R2 bases: {result['r2_bases']}")
        print(f"Total bases: {result['total_bases']}")

        if result["valid"]:
            print("Paired-end FASTQ validation passed.")
            return 0

        print("Paired-end FASTQ validation failed.", file=sys.stderr)

        for error in result["errors"]:
            print(f"Error: {error}", file=sys.stderr)

        return 1


    if args.command == "vcf-validate":
        try:
            result = validate_vcf(args.input)

        except (OSError, UnicodeError) as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1

        print(f"Valid: {result['valid']}")
        print(f"Records: {result['record_count']}")
        print(f"Errors: {result['error_count']}")

        if not result["valid"]:
            for error in result["errors"]:
                print(f"  - {error}", file=sys.stderr)
            return 1

        print("VCF validation passed.")
        return 0


    if args.command == "vcf-stats":
        try:
            stats = calculate_vcf_stats(args.input)

            write_vcf_stats(
                stats,
                args.output,
            )

        except (OSError, ValueError, UnicodeError) as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1

        print(f"Total records: {stats['total_records']}")
        print(f"Alternate alleles: {stats['total_alt_alleles']}")
        print(f"SNVs: {stats['snvs']}")
        print(f"MNVs: {stats['mnvs']}")
        print(f"Insertions: {stats['insertions']}")
        print(f"Deletions: {stats['deletions']}")
        print(f"Transitions: {stats['transitions']}")
        print(f"Transversions: {stats['transversions']}")
        print(f"Ti/Tv ratio: {stats['titv_ratio']}")
        print(f"Mean QUAL: {stats['mean_quality']}")
        print(f"Output: {args.output}")

        return 0

    
    if args.command == "vcf-filter":
        try:
            result = filter_vcf(
                input_path=args.input,
                output_path=args.output,
                min_qual=args.min_qual,
                pass_only=args.pass_only,
                chrom=args.chrom,
                min_depth=args.min_depth,
                variant_type=args.variant_type,
            )

        except (OSError, ValueError, UnicodeError) as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1

        print(f"Total records: {result['total_records']}")
        print(f"Kept records: {result['kept_records']}")
        print(f"Filtered records: {result['filtered_records']}")
        print(f"Output: {result['output_path']}")

        return 0

    
    if args.command == "vcf-extract":
        try:
            result = extract_vcf(
                input_path=args.input,
                output_path=args.output,
                chrom=args.chrom,
                start=args.start,
                end=args.end,
                variant_id=args.variant_id,
                ids_file=args.ids_file,
            )

        except (OSError, ValueError, UnicodeError) as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1

        print(f"Total records: {result['total_records']}")
        print(f"Extracted records: {result['extracted_records']}")
        print(f"Excluded records: {result['excluded_records']}")
        print(f"Output: {result['output_path']}")

        return 0



    
    return 1