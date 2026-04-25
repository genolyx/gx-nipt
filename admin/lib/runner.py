"""Subprocess helpers — run scripts and stream output to Streamlit."""
import subprocess
import sys
from pathlib import Path
from typing import Generator, List


def stream_run(cmd: List[str], cwd: Path = None) -> Generator[str, None, int]:
    """
    Yield stdout+stderr lines from a subprocess.
    Usage in Streamlit:
        placeholder = st.empty()
        log = ""
        for line in stream_run(cmd):
            log += line
            placeholder.code(log)
    """
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        cwd=str(cwd) if cwd else None,
    )
    for line in proc.stdout:
        yield line
    proc.wait()
    return proc.returncode


def run_and_collect(cmd: List[str], cwd: Path = None) -> tuple[int, str]:
    """Run a command and return (returncode, combined_output)."""
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=str(cwd) if cwd else None,
    )
    return result.returncode, result.stdout + result.stderr


def python_cmd() -> str:
    """Return the current Python executable."""
    return sys.executable
