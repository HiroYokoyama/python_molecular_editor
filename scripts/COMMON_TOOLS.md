# Shared maintenance tools

Do not edit copied scripts independently. Their tests run in the canonical repositories, not in every consumer CI matrix.

Zenodo, release and materialization helpers: [moleditpy-common-tools](https://github.com/HiroYokoyama/moleditpy-common-tools), version 1.0.0. `common-tools.lock.json` records the tested commit and LF-normalized SHA-256 for each copied script. Repository-specific JSON configuration stays beside these copies. Update from a tested checkout using `python tools/sync_consumers.py --repo python_molecular_editor --write`.
