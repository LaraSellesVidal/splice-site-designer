"""Intron analysis via Infernal (cmscan + cmalign).

Runs external Infernal tools to classify a group I intron and extract
its structural annotation, then delegates to structural_parser to
identify the P1/P1ex elements.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from .models import IntronAnalysis, P1exStemInfo
from .structural_parser import (
    StructuralElements,
    parse_cmscan_tblout,
    parse_stockholm_alignment,
)


def _find_cm_database() -> Path:
    """Locate the bundled covariance model database directory.

    Searches in order:
    1. ``GROUP_I_CM_DIR`` environment variable
    2. ``data/subtype_CMs/`` relative to this package
    3. Current working directory ``subtype_CMs/``
    """
    env_dir = os.environ.get("GROUP_I_CM_DIR")
    if env_dir:
        p = Path(env_dir)
        if p.is_dir():
            return p

    pkg_dir = Path(__file__).parent / "data" / "subtype_CMs"
    if pkg_dir.is_dir() and any(pkg_dir.glob("*.cm")):
        return pkg_dir

    cwd_dir = Path.cwd() / "subtype_CMs"
    if cwd_dir.is_dir():
        return cwd_dir

    raise FileNotFoundError(
        "Cannot find covariance model directory.  Set the GROUP_I_CM_DIR "
        "environment variable or place .cm files in data/subtype_CMs/."
    )


def _check_infernal() -> str:
    """Verify that Infernal is installed and return the path to cmscan."""
    cmscan = shutil.which("cmscan")
    if cmscan is None:
        raise RuntimeError(
            "Infernal 'cmscan' not found on PATH.  "
            "Install Infernal ≥ 1.1.4 (http://eddylab.org/infernal/)."
        )
    return cmscan


def _check_cmalign() -> str:
    """Verify that cmalign is available and return its path."""
    cmalign = shutil.which("cmalign")
    if cmalign is None:
        raise RuntimeError(
            "Infernal 'cmalign' not found on PATH.  "
            "Install Infernal ≥ 1.1.4 (http://eddylab.org/infernal/)."
        )
    return cmalign


def _write_temp_fasta(sequence: str, name: str = "query") -> str:
    """Write a sequence to a temporary FASTA file and return the path."""
    fd, path = tempfile.mkstemp(suffix=".fa", prefix="intron_")
    with os.fdopen(fd, "w") as f:
        f.write(f">{name}\n{sequence}\n")
    return path


def _run_cmscan(
    cm_dir: Path,
    fasta_path: str,
    cmscan_bin: str,
) -> tuple[str, str]:
    """Run cmscan against all .cm files in the directory.

    Returns (tblout_text, stdout_text).
    """
    # Concatenate or press CMs if needed
    cm_files = sorted(cm_dir.glob("*.cm"))
    if not cm_files:
        raise FileNotFoundError(f"No .cm files found in {cm_dir}")

    # Use a temporary tblout file
    fd, tblout_path = tempfile.mkstemp(suffix=".tblout", prefix="cmscan_")
    os.close(fd)

    try:
        # Run cmscan for each CM file individually and collect results
        all_tblout_lines: list[str] = []
        best_stdout = ""

        for cm_file in cm_files:
            # Check if the CM has been pressed (cmpress creates .i1m etc.)
            pressed = cm_file.with_suffix(".cm.i1m").exists()
            if not pressed:
                # Try to press the CM
                cmpress = shutil.which("cmpress")
                if cmpress:
                    subprocess.run(
                        [cmpress, str(cm_file)],
                        capture_output=True,
                        timeout=120,
                    )

            result = subprocess.run(
                [
                    cmscan_bin,
                    "--tblout", tblout_path,
                    "--noali",
                    "--notextw",
                    str(cm_file),
                    fasta_path,
                ],
                capture_output=True,
                text=True,
                timeout=300,
            )

            if result.returncode == 0:
                with open(tblout_path) as f:
                    for line in f:
                        if not line.startswith("#"):
                            all_tblout_lines.append(line)
                if not best_stdout:
                    best_stdout = result.stdout

        # Reconstruct a tblout-like string
        tblout_text = "\n".join(all_tblout_lines)
        return tblout_text, best_stdout

    finally:
        if os.path.exists(tblout_path):
            os.unlink(tblout_path)


def _run_cmalign(
    cm_path: Path,
    fasta_path: str,
    cmalign_bin: str,
) -> str:
    """Run cmalign and return the Stockholm alignment text."""
    result = subprocess.run(
        [
            cmalign_bin,
            "--notrunc",
            "-g",  # global alignment
            str(cm_path),
            fasta_path,
        ],
        capture_output=True,
        text=True,
        timeout=300,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"cmalign failed (exit {result.returncode}): {result.stderr}"
        )

    return result.stdout


def analyse_intron(intron_seq: str) -> IntronAnalysis:
    """Analyse a putative group I intron sequence using Infernal.

    1. Runs ``cmscan`` to classify the intron subtype.
    2. Runs ``cmalign`` against the best-matching CM to get structural annotation.
    3. Parses the alignment to identify P1, P1ex, and other elements.

    Parameters
    ----------
    intron_seq : str
        The intron DNA sequence (uppercase, no gaps).

    Returns
    -------
    IntronAnalysis
        Structural analysis results.
    """
    intron_seq = intron_seq.upper()

    cmscan_bin = _check_infernal()
    cmalign_bin = _check_cmalign()
    cm_dir = _find_cm_database()

    fasta_path = _write_temp_fasta(intron_seq)

    try:
        # Step 1: classify with cmscan
        tblout_text, _ = _run_cmscan(cm_dir, fasta_path, cmscan_bin)
        hits = parse_cmscan_tblout(tblout_text)

        if not hits:
            raise ValueError(
                "No group I intron models matched the input sequence.  "
                "Verify that the sequence is a group I intron."
            )

        best = hits[0]
        subtype = best["target_name"]
        score = best["score"]
        e_value = best["e_value"]

        # Step 2: align to best CM
        cm_path = cm_dir / f"{subtype}.cm"
        if not cm_path.exists():
            # Try to find the CM file with a different naming convention
            candidates = list(cm_dir.glob(f"*{subtype}*"))
            if candidates:
                cm_path = candidates[0]
            else:
                cm_path = sorted(cm_dir.glob("*.cm"))[0]

        stockholm = _run_cmalign(cm_path, fasta_path, cmalign_bin)

        # Step 3: parse structural annotation
        elements: StructuralElements = parse_stockholm_alignment(
            stockholm, intron_seq
        )

        # Extract native sequences
        native_igs = "".join(
            intron_seq[p] for p in elements.igs_positions
        )
        native_p1ex_p10 = "".join(
            intron_seq[p] for p in elements.p1ex_p10_positions
        )

        stem_5p_native = "".join(
            intron_seq[p] for p in elements.p1ex_stem_5prime_positions
        )
        stem_3p_native = "".join(
            intron_seq[p] for p in elements.p1ex_stem_3prime_positions
        )

        stem = P1exStemInfo(
            five_prime_positions=elements.p1ex_stem_5prime_positions,
            three_prime_positions=elements.p1ex_stem_3prime_positions,
            five_prime_native=stem_5p_native,
            three_prime_native=stem_3p_native,
            stem_length=min(
                len(elements.p1ex_stem_5prime_positions),
                len(elements.p1ex_stem_3prime_positions),
            ),
        )

        return IntronAnalysis(
            is_group_i=True,
            subtype=subtype,
            score=score,
            e_value=e_value,
            igs_positions=elements.igs_positions,
            p1ex_p10_positions=elements.p1ex_p10_positions,
            p1ex_stem=stem,
            wobble_g_position=elements.wobble_g_position,
            omega_g_position=elements.omega_g_position,
            native_igs=native_igs,
            native_p1ex_p10=native_p1ex_p10,
            native_p1_length=len(elements.igs_positions),
            native_p1ex_p10_length=len(elements.p1ex_p10_positions),
            intron_length=len(intron_seq),
        )

    finally:
        if os.path.exists(fasta_path):
            os.unlink(fasta_path)
