from __future__ import annotations

import argparse
import fnmatch
import shutil
from pathlib import Path


EXCLUDED_DIRS = {
    ".git",
    "__pycache__",
    ".pytest_cache",
    ".mypy_cache",
    ".ruff_cache",
    ".venv",
    "venv",
}

EXCLUDED_FILES = {
    ".coverage",
    "git_push.log",
}

EXCLUDED_PATTERNS = {
    "*.pyc",
    "*.pyo",
    "*.joblib",
}


def should_exclude(path: Path, root: Path, output_dir: Path) -> bool:
    if path == output_dir or output_dir in path.parents:
        return True
    relative = path.relative_to(root)
    if any(part in EXCLUDED_DIRS or part.startswith("release_bundle") for part in relative.parts):
        return True
    if path.name in EXCLUDED_FILES:
        return True
    return any(fnmatch.fnmatch(path.name, pattern) for pattern in EXCLUDED_PATTERNS)


def build_bundle(root: Path, output_dir: Path) -> tuple[int, int]:
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    copied = 0
    skipped = 0

    for path in root.rglob("*"):
        if path == output_dir:
            continue
        if should_exclude(path, root, output_dir):
            skipped += 1
            continue

        relative = path.relative_to(root)
        destination = output_dir / relative

        if path.is_dir():
            destination.mkdir(parents=True, exist_ok=True)
            continue

        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination)
        copied += 1

    return copied, skipped


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare a clean release bundle for Kaggle/HF publishing.")
    parser.add_argument(
        "--output-dir",
        default="release_bundle",
        help="Directory to write the cleaned export into.",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    output_dir = (root / args.output_dir).resolve()

    copied, skipped = build_bundle(root, output_dir)
    print(f"Prepared clean bundle at: {output_dir}")
    print(f"Copied files: {copied}")
    print(f"Skipped entries: {skipped}")


if __name__ == "__main__":
    main()
