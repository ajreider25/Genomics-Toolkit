
from itertools import zip_longest
from pathlib import Path

from .fastq_stats import read_fastq


def normalize_read_id(read_id: str) -> tuple[str, int | None]:
    """
    Return the fragment ID and optional mate number.

    Supports:
      read1/1
      read1/2
      read1 1:N:0:ATCG
      read1 2:N:0:ATCG
      read1
    """
    parts = read_id.split()
    identifier = parts[0]

    if identifier.endswith("/1"):
        return identifier[:-2], 1

    if identifier.endswith("/2"):
        return identifier[:-2], 2

    if len(parts) > 1:
        mate_field = parts[1].split(":", 1)[0]
        if mate_field in ("1", "2"):
            return identifier, int(mate_field)

    return identifier, None


def read_fastq_pairs(r1_path, r2_path):
    """
    Yield synchronized FASTQ record pairs.

    Raises ValueError when:
      - file lengths differ
      - fragment identifiers differ
      - explicit mate labels are incorrect
    """
    r1_path = Path(r1_path)
    r2_path = Path(r2_path)

    if r1_path.resolve() == r2_path.resolve():
        raise ValueError("R1 and R2 must be different files.")

    if not r1_path.is_file():
        raise FileNotFoundError(f"R1 file not found: {r1_path}")

    if not r2_path.is_file():
        raise FileNotFoundError(f"R2 file not found: {r2_path}")

    for pair_number, (r1, r2) in enumerate(
        zip_longest(read_fastq(r1_path), read_fastq(r2_path)),
        start=1,
    ):
        if r1 is None or r2 is None:
            raise ValueError(
                f"R1 and R2 have different read counts "
                f"(first unmatched pair: {pair_number})."
            )

        r1_id, r1_mate = normalize_read_id(r1[0])
        r2_id, r2_mate = normalize_read_id(r2[0])

        if r1_id != r2_id:
            raise ValueError(
                f"Pair {pair_number}: mismatched read IDs: "
                f"{r1[0]!r} vs {r2[0]!r}"
            )

        if r1_mate is not None and r1_mate != 1:
            raise ValueError(
                f"Pair {pair_number}: R1 has incorrect mate label."
            )

        if r2_mate is not None and r2_mate != 2:
            raise ValueError(
                f"Pair {pair_number}: R2 has incorrect mate label."
            )

        yield r1, r2


def validate_fastq_pairs(r1_path, r2_path) -> dict:
    """Validate paired FASTQ files and return a summary."""
    pair_count = 0
    r1_bases = 0
    r2_bases = 0
    errors = []

    try:
        for r1, r2 in read_fastq_pairs(r1_path, r2_path):
            pair_count += 1
            r1_bases += len(r1[1])
            r2_bases += len(r2[1])

    except (OSError, ValueError, UnicodeError) as exc:
        errors.append(str(exc))

    return {
        "valid": not errors,
        "pair_count": pair_count,
        "r1_bases": r1_bases,
        "r2_bases": r2_bases,
        "total_bases": r1_bases + r2_bases,
        "error_count": len(errors),
        "errors": errors,
    }
