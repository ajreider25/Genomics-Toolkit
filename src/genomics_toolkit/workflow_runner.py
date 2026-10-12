import subprocess
import sys
from pathlib import Path

def run_workflow(
    cores: int = 1,
    snakefile: str | Path = "workflow/Snakefile",
    input_path: str | Path | None = None,
    output_path: str | Path | None = None
) -> int:
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
        "--snakefile",
        str(snakefile),
        "--cores",
        str(cores)
    ]

    workflow_config = {}

    if input_path is not None:
        workflow_config["input"] = str(
            Path(input_path).resolve()
        )

    if output_path is not None:
        workflow_config["output"] = str(
            Path(output_path).resolve()
        )

    if workflow_config:
        command.append("--config")

        command.extend(
            f"{key}={value}"
            for key, value in workflow_config.items()
        )

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