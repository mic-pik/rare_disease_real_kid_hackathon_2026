#!/usr/bin/env python3
"""Download an explicitly authorised gated dataset; credentials stay outside code."""
import argparse
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--dataset', required=True)
    ap.add_argument('--revision', required=True, help='Exact dataset commit')
    ap.add_argument('--out', required=True, type=Path)
    args = ap.parse_args()
    from huggingface_hub import snapshot_download
    snapshot_download(args.dataset, repo_type='dataset', revision=args.revision, local_dir=args.out)


if __name__ == '__main__': main()
