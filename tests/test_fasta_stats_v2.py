
import gzip

from genomics_toolkit.operations.fasta_stats import (
    calculate_stats,
    fasta_stats,
)


def test_nucleotide_composition():
    result = calculate_stats("AACCGGTTNN")

    assert result["length"] == 10
    assert result["a_count"] == 2
    assert result["c_count"] == 2
    assert result["g_count"] == 2
    assert result["t_count"] == 2
    assert result["n_count"] == 2
    assert result["gc_percent"] == 50.0
    assert result["at_percent"] == 50.0


def test_lowercase_nucleotide_composition():
    result = calculate_stats("aaacgtn")

    assert result["a_count"] == 3
    assert result["c_count"] == 1
    assert result["g_count"] == 1
    assert result["t_count"] == 1
    assert result["n_count"] == 1
    assert result["gc_percent"] == 33.33
    assert result["at_percent"] == 66.67


def test_ambiguous_bases_excluded_from_percentages():
    result = calculate_stats("GCNNRY")

    assert result["length"] == 6
    assert result["gc_percent"] == 100.0
    assert result["at_percent"] == 0.0
    assert result["n_count"] == 2


def test_only_ambiguous_bases():
    result = calculate_stats("NNRY")

    assert result["gc_percent"] == 0.0
    assert result["at_percent"] == 0.0


def test_gzipped_fasta(tmp_path):
    fasta_file = tmp_path / "sample.fasta.gz"

    with gzip.open(fasta_file, "wt", encoding="utf-8") as handle:
        handle.write(
            ">seq1\n"
            "AACCGGTT\n"
            ">seq2\n"
            "GGCCNN\n"
        )

    results = list(fasta_stats(fasta_file))

    assert len(results) == 2

    assert results[0]["sequence_id"] == "seq1"
    assert results[0]["length"] == 8
    assert results[0]["gc_percent"] == 50.0
    assert results[0]["at_percent"] == 50.0

    assert results[1]["sequence_id"] == "seq2"
    assert results[1]["length"] == 6
    assert results[1]["gc_percent"] == 100.0
    assert results[1]["n_count"] == 2


def test_empty_gzipped_fasta(tmp_path):
    fasta_file = tmp_path / "empty.fasta.gz"

    with gzip.open(fasta_file, "wt", encoding="utf-8"):
        pass

    assert list(fasta_stats(fasta_file)) == []
