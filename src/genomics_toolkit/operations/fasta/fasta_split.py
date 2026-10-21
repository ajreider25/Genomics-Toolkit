import gzip
import shutil
import tempfile
from pathlib import Path

from .fasta_stats import read_fasta

def split_fasta(
    input_path: str | Path,
    output_dir: str | Path,
    records_per_file: int = 1000,
    gzip_output: bool = False
) -> dict:
    """Split a FASTA into file containing at most N records."""
    input_path = Path(input_path)
    output_dir = Path(output_dir)

    if not input_path.is_file():
        raise FileExistsError(
            f"FASTA input file found: {input_path}"
        )

    if records_per_file < 1:
        raise ValueError("records_per_file must at least 1.")

    if output_dir.exists():
        raise FileExistsError(
            f"Output directory already exists: {output_dir}"
        )

    # Avoid placing output in a directory with the same path
    # as the input file.
    if input_path.resolve() == output_dir.resolve():
        raise ValueError(
            "Input file and output directory must differ."
        )

    output_dir.parent.mkdir(parents=True, exist_ok=True)

    # Build chuncks in a temporary directory, then publish
    # the directory only after processing succeeds.
    temporary_dir = Path(
        tempfile.mkdtemp(
            prefix=f".{output_dir.name}.",
            dir=output_dir.parent
        )
    )

    total_records = 0
    file_count = 0
    handle = None

    try:
        for sequence_id, sequence in read_fasta(input_path):
            if total_records % records_per_file == 0:
                if handle is not None:
                    handle.close()

                file_count += 1

                extension = (
                    ".fasta.gz" if gzip_output else ".fasta"
                )
                filename = f"part_{file_count:03d}{extension}"
                part_path = temporary_dir / filename

                opener = gzip.open if gzip_output else open
                handle = opener(
                    part_path,
                    "wt",
                    encoding="utf-8"
                )

            handle.write(f">{sequence_id}\n")

            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")

            total_records += 1

        if handle is not None:
            handle.close()
            handle = None

        # Publish completed output directory.
        # Rejecting existing destinations avoids overwrites.
        temporary_dir.rename(output_dir)

    finally:
        if handle is not None:
            handle.close()

        if temporary_dir.exists():
            shutil.rmtree(temporary_dir)

    return {
        "total_records": total_records,
        "file_count": file_count,
        "records_per_file": records_per_file,
        "output_dir": str(output_dir)
    }