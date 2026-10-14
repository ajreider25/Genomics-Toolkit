import csv
import gzip
from collections import Counter
from pathlib import Path
from typing import Iterator

def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield (sequence_id, sequence) pairs from a FASTA file."""

    path = Path(path)

    opener = gzip.open if path.suffix.lower() == ".gz" else open

    sequence_id = None
    sequence_parts = []

    with opener(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if sequence_id is not None:
                    yield sequence_id, "".join(sequence_parts)

                sequence_id = line[1:].split()[0] if line[1:].strip() else ""
                if not sequence_id:
                    raise ValueError(
                        f"Invalid FASTA header at line {line_number}"
                    )
                
                sequence_parts = []
            else:
                if sequence_id is None:
                    raise ValueError(
                        f"Invalid FASTA: sequence before header at line {line_number}"
                    )
                sequence_parts.append(line.upper())

    if sequence_id is not None:
        yield sequence_id, "".join(sequence_parts)

def calculate_stats(sequence: str) -> dict:
    """Calculate basic stats for a DNA sequence"""
    sequence = sequence.upper()

    a_count = sequence.count("A")
    c_count = sequence.count("C")
    g_count = sequence.count("G")
    t_count = sequence.count("T")
    n_count = sequence.count("N")

    canonical_bases = a_count + c_count + g_count + t_count
    gc_count = g_count + c_count
    at_count = a_count + t_count

    gc_percent = (
        round(gc_count / canonical_bases * 100, 2)
        if canonical_bases
        else 0.0
    )

    at_percent = (
        round(at_count / canonical_bases * 100, 2)
        if canonical_bases
        else 0.0
    )

    return {
        "length": len(sequence),
        "gc_percent": gc_percent,
        "at_percent": at_percent,
        "n_count": n_count,
        "a_count": a_count,
        "c_count": c_count,
        "g_count": g_count,
        "t_count": t_count
    }

def fasta_stats(path: str | Path) -> Iterator[dict]:
    """Generate stats for every FASTA record."""

    for sequence_id, sequence in read_fasta(path):
        yield {
            "sequence_id": sequence_id,
            **calculate_stats(sequence)
        }

def write_fasta_stats(
        input_path: str | Path,
        output_path: str | Path
) -> None:
    """Write FASTA sequence stats to a CSV file"""

    input_path = Path(input_path)
    output_path = Path(output_path)

    if not input_path.is_file():
        raise FileNotFoundError(
            f"FASTA input file not found: {input_path}"
        )
    
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    fieldnames = [
        "sequence_id",
        "length",
        "gc_percent",
        "at_percent",
        "n_count",
        "a_count",
        "c_count",
        "g_count",
        "t_count"
    ]

    with output_path.open(
        "w", newline="", encoding="utf-8"
    ) as file:
        writer = csv.DictWriter(
            file,
            fieldnames=fieldnames
        )

        writer.writeheader()

        for record in fasta_stats(input_path):
            writer.writerow(record)