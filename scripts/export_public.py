#!/usr/bin/env python3
"""Export the allowlisted public projection into an independent Git repo."""

from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
ALLOWLIST = ROOT / "PUBLIC_FILES.txt"
MANIFEST_NAME = ".public-export-manifest.json"
CHECKSUM_NAME = "CHECKSUMS.sha256"
FORBIDDEN_PARTS = {
    ".git", ".env", "private", "secrets", "tmp", "temp",
}
FORBIDDEN_FRAGMENTS = {
    "secret", "credential", "token", "cookie", "session", "private",
}


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def read_allowlist() -> list[Path]:
    paths: list[Path] = []
    for raw in ALLOWLIST.read_text().splitlines():
        entry = raw.strip()
        if not entry or entry.startswith("#"):
            continue
        relative = Path(entry)
        if relative.is_absolute() or ".." in relative.parts:
            raise SystemExit(f"unsafe allowlist path: {entry}")
        lowered = {part.lower() for part in relative.parts}
        contains_forbidden_fragment = any(
            fragment in part
            for part in lowered
            for fragment in FORBIDDEN_FRAGMENTS
        )
        contains_env_file = any(
            part == ".env" or part.startswith(".env.") for part in lowered
        )
        if lowered & FORBIDDEN_PARTS or contains_forbidden_fragment or contains_env_file:
            raise SystemExit(f"forbidden allowlist path: {entry}")
        source = ROOT / relative
        if not source.is_file() or source.is_symlink():
            raise SystemExit(f"allowlisted path is not a regular file: {entry}")
        paths.append(relative)
    if len(paths) != len(set(paths)):
        raise SystemExit("duplicate path in PUBLIC_FILES.txt")
    return sorted(paths, key=lambda path: path.as_posix())


def load_previous(target: Path) -> dict[str, str] | None:
    manifest = target / MANIFEST_NAME
    if not manifest.exists():
        return None
    payload = json.loads(manifest.read_text())
    return dict(payload.get("files", {}))


def ensure_safe_target(target: Path) -> None:
    resolved = target.resolve()
    if resolved in {Path("/").resolve(), Path.home().resolve(), ROOT.resolve()}:
        raise SystemExit(f"refusing unsafe export target: {resolved}")
    if not (resolved / ".git").exists():
        raise SystemExit("export target must be an independent initialized Git repository")


def assert_managed_file_unchanged(target: Path, relative: Path, old_hash: str) -> None:
    destination = target / relative
    if destination.exists() and (not destination.is_file() or digest(destination) != old_hash):
        raise SystemExit(f"refusing to overwrite modified public file: {relative}")


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: export_public.py /absolute/path/to/public-repository")
    target = Path(sys.argv[1]).expanduser()
    if not target.is_absolute():
        raise SystemExit("export target must be an absolute path")
    ensure_safe_target(target)

    allowlisted = read_allowlist()
    previous = load_previous(target)
    current_names = {path.as_posix() for path in allowlisted}

    if previous is None:
        collisions = [relative for relative in allowlisted
                      if (target / relative).exists()]
        if collisions:
            rendered = ", ".join(path.as_posix() for path in collisions)
            raise SystemExit(
                "first export requires an empty public working tree; "
                f"existing allowlisted paths: {rendered}"
            )
        previous = {}
    else:
        for relative_text, old_hash in previous.items():
            assert_managed_file_unchanged(target, Path(relative_text), old_hash)

    for relative in allowlisted:
        destination = target / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / relative, destination)

    for relative_text, old_hash in previous.items():
        if relative_text not in current_names:
            destination = target / relative_text
            if destination.exists() and digest(destination) == old_hash:
                destination.unlink()

    hashes = {relative.as_posix(): digest(target / relative)
              for relative in allowlisted}
    checksum_lines = [f"{hashes[name]}  {name}" for name in sorted(hashes)]
    (target / CHECKSUM_NAME).write_text("\n".join(checksum_lines) + "\n")
    manifest = {
        "schema": "m23.public-export.v1",
        "source_policy": "allowlist",
        "files": hashes,
        "checksum_file": CHECKSUM_NAME,
    }
    (target / MANIFEST_NAME).write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"exported {len(hashes)} managed files to {target}")
    print(f"review {MANIFEST_NAME}, {CHECKSUM_NAME}, git status, and git diff")


if __name__ == "__main__":
    main()
