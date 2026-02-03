#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import yaml

AUTOGEN_SENTINEL = "# AUTO-GENERATED FILE. DO NOT EDIT BY HAND."

MANUAL_BEGIN = "# === MANUAL BEGIN ==="
MANUAL_END   = "# === MANUAL END ==="


def _extract_manual_block(text: str) -> str | None:
    """
    Return the manual block including markers, or None if not present/invalid.
    """
    i = text.find(MANUAL_BEGIN)
    j = text.find(MANUAL_END)
    if i == -1 or j == -1 or j < i:
        return None
    j_end = j + len(MANUAL_END)
    return text[i:j_end].rstrip() + "\n"


def _inject_manual_block(text: str, block: str) -> str:
    """
    Replace existing manual block if present, otherwise append at end.
    """
    existing = _extract_manual_block(text)
    if existing is not None:
        return text.replace(existing, block)
    # append with a blank line before if needed
    if not text.endswith("\n"):
        text += "\n"
    if not text.endswith("\n\n"):
        text += "\n"
    return text + block


def load_yaml(p: Path) -> dict:
    return yaml.safe_load(p.read_text(encoding="utf-8")) or {}


def is_autogen_file(p: Path) -> bool:
    if not p.exists():
        return True
    head = p.read_text(encoding="utf-8", errors="ignore").splitlines()[:25]
    return any(AUTOGEN_SENTINEL in line for line in head)


def safe_copy(src: Path, dst: Path, *, enforce_autogen: bool) -> None:
    if enforce_autogen and dst.exists() and not is_autogen_file(dst):
        raise SystemExit(f"Refusing to overwrite non-autogen file: {dst}")

    dst.parent.mkdir(parents=True, exist_ok=True)

    # --- NEW: preserve manual block from existing dst (if any) ---
    manual_block = None
    if dst.exists() and is_autogen_file(dst):
        dst_text = dst.read_text(encoding="utf-8", errors="ignore")
        manual_block = _extract_manual_block(dst_text)

    src_text = src.read_text(encoding="utf-8", errors="strict")
    if manual_block is not None:
        src_text = _inject_manual_block(src_text, manual_block)

    tmp = dst.with_suffix(dst.suffix + ".tmp")
    tmp.write_text(src_text, encoding="utf-8")
    tmp.replace(dst)


def copy_package_folder(
    *,
    src_pkg_dir: Path,
    dst_pkg_dir: Path,
    target_filename: str,
    copy_init: bool,
    enforce_autogen: bool,
    repo_root: Path,
    dry_run: bool,
) -> tuple[bool, bool]:
    """
    Copy <src_pkg_dir>/<target_filename> and optionally <src_pkg_dir>/__init__.py
    into <dst_pkg_dir>/...
    Returns (copied_models, copied_init).
    """
    src_models = src_pkg_dir / target_filename
    if not src_models.exists():
        raise SystemExit(f"Missing {target_filename} in source folder: {src_pkg_dir}")

    dst_models = dst_pkg_dir / target_filename
    if dry_run:
        print(f"[copy] {src_models.relative_to(repo_root)} -> {dst_models.relative_to(repo_root)}")
    else:
        safe_copy(src_models, dst_models, enforce_autogen=enforce_autogen)
        print(f"[copy] {src_models.relative_to(repo_root)} -> {dst_models.relative_to(repo_root)}")

    copied_init = False
    if copy_init:
        src_init = src_pkg_dir / "__init__.py"
        if src_init.exists():
            dst_init = dst_pkg_dir / "__init__.py"
            if dry_run:
                print(f"[copy] {src_init.relative_to(repo_root)} -> {dst_init.relative_to(repo_root)}")
            else:
                safe_copy(src_init, dst_init, enforce_autogen=enforce_autogen)
                print(f"[copy] {src_init.relative_to(repo_root)} -> {dst_init.relative_to(repo_root)}")
            copied_init = True

    return True, copied_init


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, required=True, help="sync_offbeat_autogen_map.yaml")
    ap.add_argument("--repo-root", type=Path, default=Path("."), help="Repo root (default: .)")
    ap.add_argument("--dry-run", action="store_true", help="Print actions without writing")
    ap.add_argument("--allow-missing", action="store_true", help="Skip missing mapped source folders")
    args = ap.parse_args()

    repo_root = args.repo_root.resolve()
    cfg = load_yaml(args.config)

    source_dir = (repo_root / cfg["source_dir"]).resolve()
    dest_root = (repo_root / cfg["dest_root"]).resolve()

    target_filename: str = str(cfg.get("target_filename", "models.py")).strip() or "models.py"
    copy_init: bool = bool(cfg.get("copy_init", True))
    enforce_autogen: bool = bool(cfg.get("enforce_autogen", False))

    folder_map: dict[str, str] = cfg.get("folder_map", {}) or {}
    if not isinstance(folder_map, dict) or not folder_map:
        raise SystemExit("Config must define non-empty 'folder_map' mapping source folders -> dest folders")

    if not source_dir.exists():
        raise SystemExit(f"source_dir not found: {source_dir}")
    if not dest_root.exists():
        raise SystemExit(f"dest_root not found: {dest_root}")

    n_models = 0
    n_inits = 0

    # Copy only what is explicitly mapped
    for src_rel_str, dst_rel_str in folder_map.items():
        src_rel = Path(str(src_rel_str))
        dst_rel = Path(str(dst_rel_str))

        src_pkg_dir = (source_dir / src_rel).resolve()
        dst_pkg_dir = (dest_root / dst_rel).resolve()

        if not src_pkg_dir.exists():
            msg = f"Mapped source folder not found: {src_pkg_dir}  (map key: {src_rel_str!r})"
            if args.allow_missing:
                print(f"[skip] {msg}")
                continue
            raise SystemExit(msg)

        copied_models, copied_init = copy_package_folder(
            src_pkg_dir=src_pkg_dir,
            dst_pkg_dir=dst_pkg_dir,
            target_filename=target_filename,
            copy_init=copy_init,
            enforce_autogen=enforce_autogen,
            repo_root=repo_root,
            dry_run=args.dry_run,
        )
        if copied_models:
            n_models += 1
        if copied_init:
            n_inits += 1

    print(f"\nDone. Copied {n_models} '{target_filename}' files and {n_inits} __init__.py files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
