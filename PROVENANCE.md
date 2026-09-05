# Provenance and scope

This publication copy assembles three successive locally developed benchmark
packages, labeled v0.1, v0.2 and v0.3. These are research-checkpoint labels, not
software stability guarantees or peer-review certifications.

The equations, computational algorithms, fixed numerical seeds, raw numerical
result values, error analyses, and third-party research citations are retained.
Source module filenames and their one internal import were renamed. Headings,
unpublished-project cross-references, and local time metadata were edited for a
neutral presentation. The original private archive checksums are not published;
SHA256SUMS.txt covers the contents of this publication copy instead.

All four Python sources were compared to their predecessors as abstract syntax
trees, excluding module docstrings and the renamed internal import. Their
executable syntax was otherwise unchanged. Numeric seeds may look like dates;
they remain arbitrary reproducibility inputs, not certified run timestamps.

Original high-repetition numerical results remain in each stage's results.json
and CSV files. The new packaging verification is separately labeled and does not
claim to rerun every original high-repetition simulation. Normal reproduction
writes to outputs/, not over the original records.

AI assistance was used in developing the derivations, simulations and text.
The benchmark tests and packaging checks are not an independent replication,
a hardware experiment, or peer review. No claim of institutional affiliation
or external endorsement is made.

Removing explicit identity-related metadata does not make code or findings
untraceable. The preserved scientific content, wording, results and seeds can be
matched to other copies. This package offers no anonymity guarantee.
