import subprocess
import sys
from pathlib import Path

def run_workflow(cores: int = 1, snakefile: str | Path = "workflow/Snakefile") -> int:
    """Execute a Snakemake workflow."""
    snakefile = Path(snakefile).resolve()

    if not snakefile.is_file():
        print(
            f"Error: Snakefile not found: {snakefile}",
            file=sys.stderr
        )
        return 1

    command = [
        sys.executable,
        "-m",
        "snakemake",
        "--snakemake",
        str(snakefile),
        "--cores",
        str(cores)
    ]

    try:
        result = subprocess.run(
            command,
            check=False
        )
    except OSError as exc:
        print(
            f"Failed to start Snakemake: {exc}",
            file=sys.stderr
        )
        return 1

    return result.returncode