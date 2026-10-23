"""Reuse the same bootstrap runner with longer blocks and frozen main fits."""
import argparse
import json
import os
from pathlib import Path
from dispersion_run import main
from io_utils import digest, dump


def run(source, reference, output, session):
    ref = (reference/session/'main').resolve()
    dest = output/session; dest.mkdir(parents=True, exist_ok=True)
    signature = {p.name:digest(p) for p in sorted(ref.glob('*.json'))}
    stamp = dest/'reference_signature.json'
    if stamp.exists() and json.loads(stamp.read_text()) != signature:
        raise ValueError('Reference fits changed')
    dump(stamp, signature)
    link = dest/'main'
    if link.exists():
        if link.resolve() != ref:
            raise ValueError('Wrong reference directory')
    else:
        link.symlink_to(os.path.relpath(ref, dest.resolve()), target_is_directory=True)
    main(source, output, session, Path(__file__).with_name('dispersion_block10.toml'), 'bootstrap')


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--reference',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--session',choices=['590','n0','20w'],required=True)
    a=p.parse_args();run(a.source,a.reference,a.output,a.session)
