import os
import gzip
import tempfile
from pathlib import Path

from .fasta_stats import read_fasta

COMPLEMENT_TABLE = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)

VALID_DNA_BASES = set("ACGTRYSWKMBDHVN")

def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    invalid =set(sequence.upper()) - VALID_DNA_BASES

    if invalid:
        raise ValueError(
            "Invalid DNA nucleotide(s): "
            + ", ".join(sorted(invalid))
        )

    return sequence.translate(COMPLEMENT_TABLE)[::-1]

def write_reverse_complement(
    input_path: str | Path,
    output_path: str | Path
) -> dict:
    """Write reverse-complemented FASTA records atomically."""
    input_path = Path(input_path)
    output_path = Path(output_path)

    if not input_path.is_file():
        raise FileNotFoundError(
            f"FASTA input file not found: {input_path}"
        )

    if input_path.resolve() == output_path.resolve():
        raise ValueError(
            "Input and output cannot be the same file."
        )
    
    output_path.parent.mkdir(parents=True, exist_ok=True)

    sequence_count = 0
    total_bases = 0
    temporary_path = None

    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=output_path.parent,
            prefix=f".{output_path.name}.",
            suffix=".tmp",
            delete=False                
        ) as temporary_file:
            temporary_path = Path(temporary_file.name)

        opener = gzip.open if output_path.suffix == ".gz" else open

        with opener(
            temporary_path,
            "wt",
            encoding="utf-8"
        ) as handle:
            for sequence_id, sequence in read_fasta(input_path):
                revcomp = reverse_complement(sequence)

                handle.write(f">{sequence_id}\n")

                for start in range(0, len(revcomp), 80):
                    handle.write(revcomp[start:start + 80] + "\n")

                sequence_count +=1
                total_bases += len(sequence)

        os.replace(temporary_path, output_path)
        temporary_path = None

    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)

    return {
        "sequence_count": sequence_count,
        "total_bases": total_bases
    }