#!/usr/bin/env python3
"""Scan distributable files without printing potential secret values."""
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN_SUFFIXES = {'.pt','.pth','.safetensors','.npz','.npy','.parquet','.joblib','.bam','.cram','.pem','.key'}
SECRET_PATTERNS = [
    re.compile(r'\bh' + r'f_[A-Za-z0-9]{20,}\b'),
    re.compile(r'\bs' + r'k-(?:proj-|ant-)?[A-Za-z0-9_-]{25,}\b'),
    re.compile(r'\bg' + r'h[pousr]_[A-Za-z0-9]{25,}\b'),
    re.compile(r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----'),
    re.compile(r'(?i)(?:api_key|access_token|secret_key)\s*=\s*[\"\x27][A-Za-z0-9_-]{20,}')]
MACHINE_PATH = re.compile(r'/(?:data/shared|home|Users)/[A-Za-z0-9_.-]+')


def candidates():
    # Include tracked files even if ignored, but omit files deleted from the new tree.
    names = subprocess.check_output(['git','ls-files','--cached','--others','--exclude-standard','-z'], cwd=ROOT)
    return [ROOT/x.decode() for x in sorted(set(names.split(b'\0'))) if x and (ROOT/x.decode()).is_file()]


def main():
    problems = []
    files = candidates()
    for p in files:
        rel = p.relative_to(ROOT)
        if p.is_symlink(): problems.append((str(rel), 'symlink is not distributable'))
        if p.suffix in FORBIDDEN_SUFFIXES or any(x in p.name.lower() for x in ('.vcf','.fastq','.cram','.bam')):
            problems.append((str(rel), 'data or checkpoint extension'))
        if p.stat().st_size > 5_000_000: problems.append((str(rel), 'unexpected large artifact'))
        try: text = p.read_text()
        except UnicodeDecodeError:
            problems.append((str(rel), 'binary file')); continue
        for i, line in enumerate(text.splitlines(), 1):
            if MACHINE_PATH.search(line): problems.append((f'{rel}:{i}', 'machine-specific path'))
            if any(pattern.search(line) for pattern in SECRET_PATTERNS):
                problems.append((f'{rel}:{i}', 'possible credential; value suppressed'))
    for location, issue in problems: print(f'{location}: {issue}')
    print(f'Checked {len(files)} files; {len(problems)} findings. No secret values printed.')
    raise SystemExit(bool(problems))


if __name__ == '__main__': main()
