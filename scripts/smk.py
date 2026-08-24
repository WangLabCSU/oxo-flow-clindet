#!/usr/bin/env python3
"""smk.py — shim that lets verbatim upstream clindet Python scripts (which
talk to Snakemake through the global `snakemake` module) run under oxo-flow.

Usage:
  python3 scripts/smk.py --script scripts/merge_caller_vcfs.py \\
      --input vcfs f1 f2 ... --param dir DIR --output merged_vcf OUT

Flags (same convention as scripts/smk.R):
  --script PATH         Python script to run (exec'd verbatim)
  --input NAME V...     named input slot (multiple values = list)
  --output NAME V...    named output slot
  --param NAME V...     named params slot
  --wildcard NAME V...  named wildcards slot
  --threads N           threads slot

The shim inserts a module named `snakemake` into sys.modules with the same
attribute surface the upstream scripts touch (input/output/params/wildcards/
threads/log), then execs the script.
"""
import argparse
import sys
import types


class _Slot:
    def __init__(self, scalar_when_single=False):
        self._items = []
        self._scalar_when_single = scalar_when_single

    def _resolve(self, name):
        # Upstream scripts read snakemake.input.* as lists (expand(...))
        # and snakemake.output.*/params.* as scalars.
        values = [v for (n, v) in self._items if n == name]
        if not values:
            raise KeyError(name)
        flat = values[0] if len(values) == 1 else [
            item for vals in values for item in vals
        ]
        if self._scalar_when_single and len(flat) == 1:
            return flat[0]
        return flat

    def __getattr__(self, name):
        try:
            return self._resolve(name)
        except KeyError:
            raise AttributeError(
                f"snakemake slot '{name}' was not provided"
            ) from None

    # dict-style access (upstream config_freec.py uses
    # snakemake.input['ini_template'] / snakemake.input.keys()): single-file
    # slots resolve to scalars, matching real Snakemake.
    def __getitem__(self, name):
        values = [v for (n, v) in self._items if n == name]
        if not values:
            raise KeyError(name)
        flat = values[0] if len(values) == 1 else [
            item for vals in values for item in vals
        ]
        if len(flat) == 1:
            return flat[0]
        return flat

    def keys(self):
        return [n for (n, _v) in self._items]

    def __contains__(self, name):
        return any(n == name for (n, _v) in self._items)

    def __setattr__(self, name, value):
        if name in ("_items", "_scalar_when_single"):
            object.__setattr__(self, name, value)
        else:
            self._items.append((name, value))


def build_module(args):
    mod = types.ModuleType("snakemake")
    mod.input = _Slot(scalar_when_single=False)
    mod.output = _Slot(scalar_when_single=True)
    mod.params = _Slot(scalar_when_single=True)
    mod.wildcards = _Slot(scalar_when_single=True)
    mod.log = _Slot(scalar_when_single=True)
    mod.threads = 1

    script = None
    i = 0
    while i < len(args):
        flag = args[i][2:]
        if flag == "script":
            script = args[i + 1]
            i += 2
            continue
        if flag == "threads":
            mod.threads = int(args[i + 1])
            i += 2
            continue
        if flag in ("input", "output", "param", "wildcard", "log"):
            name = args[i + 1]
            vals = []
            j = i + 2
            while j < len(args) and not args[j].startswith("--"):
                vals.append(args[j])
                j += 1
            slot = getattr(mod, flag if flag != "param" else "params")
            # All the upstream scripts read input slots as lists
            # (snakemake.input.vcfs = expand(...)), so keep list semantics.
            slot._items.append((name, vals))
            i = j
            continue
        i += 1
    return mod, script


def main():
    # argparse would mangle --input values; parse manually (smk.R convention).
    mod, script = build_module(sys.argv[1:])
    if script is None:
        sys.exit("smk.py: --script PATH is required")
    sys.modules["snakemake"] = mod
    with open(script) as fh:
        code = compile(fh.read(), script, "exec")
        exec(code, {"__name__": "__main__", "__file__": script})


if __name__ == "__main__":
    main()
