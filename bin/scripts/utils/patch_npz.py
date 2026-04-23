#!/usr/bin/env python3
"""Patch old WisecondorX reference NPZ files to add trained_cutoff and is_nipt.

WisecondorX 1.2+ requires these keys which are absent in references built
with older versions (1.1.x and below).
"""
import sys, os, glob
import numpy as np

TRAINED_CUTOFF = np.float64(0.003)  # Y/total reads threshold for NIPT gender
IS_NIPT        = np.bool_(True)

def patch_npz(path: str) -> str:
    d = np.load(path, allow_pickle=True)
    keys = list(d.keys())
    if 'trained_cutoff' in keys:
        return 'skip'
    data = {}
    for k in keys:
        if k == 'arguments':
            continue   # old pickled toolNewref object – not needed by predict
        try:
            data[k] = d[k]
        except Exception:
            pass
    data['trained_cutoff'] = TRAINED_CUTOFF
    data['is_nipt']        = IS_NIPT
    # np.savez_compressed auto-appends .npz, so use a tmp path without .npz
    tmp_base = path[:-4] + '.patching'   # strip .npz, add different suffix
    np.savez_compressed(tmp_base, **data)  # writes tmp_base + '.npz'
    tmp_path = tmp_base + '.npz'
    os.replace(tmp_path, path)
    return 'patched'


if __name__ == '__main__':
    search_root = sys.argv[1] if len(sys.argv) > 1 else '.'
    files = glob.glob(os.path.join(search_root, '**', '*.npz'), recursive=True)
    print(f"Found {len(files)} NPZ files under {search_root}")
    p = s = f = 0
    for path in sorted(files):
        try:
            result = patch_npz(path)
            if result == 'patched':
                print(f"  PATCHED: {path}")
                p += 1
            else:
                s += 1
        except Exception as e:
            print(f"  FAILED:  {path} ({e})")
            f += 1
    print(f"\nDone: patched={p}  skipped={s}  failed={f}")
