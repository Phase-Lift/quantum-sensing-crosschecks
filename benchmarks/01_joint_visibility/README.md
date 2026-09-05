# Joint visibility benchmark v0.1

See CHECKPOINT.md for the GHZ construction, proper-subset blindness, exact
classification calculations, noise assumptions, and resource accounting.

```sh
python joint_visibility.py --trials 20000 --output results.json
```

Running in this directory with that output name replaces its recorded results.
Use the repository-root run_checks.py to preserve those files. This is a simulation,
not a hardware run or a sub-Planck observation.
