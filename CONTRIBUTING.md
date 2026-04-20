# Contributing to thermocas

## Development setup

```bash
git clone <repo>
cd thermocas-framework
uv venv && source .venv/bin/activate
uv pip install -e ".[dev]"
pytest
```

## Project layout

```
src/thermocas/        # framework code
tests/                # pytest suite — one test file per module
examples/             # runnable demos (synthetic_pipeline.py, v1/v2_pipeline.sh)
config/               # PAM model + cohort YAMLs
data/  results/       # gitignored — local artifacts only
```

### Reproducible artifacts vs primary data

`data/derived/scored_*.jsonl` and `data/derived/catalog_*.jsonl` are
**reproducible outputs**, not primary assets, and are intentionally
gitignored. They are derivable from the committed summary TSVs
(`data/derived/*_cohort/tumor_summary.tsv` +
`data/derived/*_cohort/normal_summary.tsv` +
`data/derived/*_cohort/probes.tsv`) plus the PAM model via:

```bash
thermocas build-catalog  ...   # regenerate the catalog JSONL
thermocas score-cohort   ...   # regenerate a scored cohort JSONL
```

These artifacts were removed from git history via `git-filter-repo` on
2026-04-20 to keep the repo small and cloneable; any reference to
pre-rewrite commit SHAs is stale. If you have a local clone from before
that rewrite, reclone rather than merge.

## Conventions

- **Stdlib only** for the core framework (no numpy / scipy / pandas / requests).
  Optional accelerators can ship in a separate extras group later.
- **Pydantic v2 models** are the contract. Every persisted record (catalog,
  scored cohort, pan-cancer atlas, benchmark result) is one Pydantic model
  per JSONL line.
- **Validators fail fast.** New invariants belong in `models.py` validators,
  not in downstream defensive checks.
- **Streaming over batching.** New helpers should yield rather than return
  full lists where genome scale is plausible.
- **Score components are auditable.** Don't return opaque scalars; return
  decomposed Pydantic models so reviewers can see why a candidate ranked.

## Running tests

```bash
pytest                    # full suite, ~1s
pytest tests/test_grna.py # one module
pytest -k probabilistic   # by name
ruff check src tests      # linting
```

The full suite must pass before sending a PR. CI runs the same commands.

## Adding a new module

1. Add a Pydantic model to `models.py` for any persisted record.
2. Add the module under `src/thermocas/` with a one-paragraph header
   explaining what problem it solves and how it composes with existing layers.
3. Add tests under `tests/` covering: happy path, every public boundary,
   one explicit edge case per invariant.
4. If the module ships a CLI surface, add a subparser in `cli.py` and a
   CLI test in `tests/test_cli.py` that runs the full subcommand against
   a tmp_path fixture.
5. Update `CHANGELOG.md` under the **Unreleased** heading.
6. Update the V*-status table in `README.md`.

## Commit messages

Conventional-commits style. Examples:

```
feat(probabilistic): add Beta-distribution CDF
fix(pam_model): reject mixed-width regexes at load time
docs(README): describe gdc-fetch subcommand
test(grna): cover hairpin-rich periodic spacers
```

## Review checklist

Before sending a PR, verify:

- [ ] `pytest` is green
- [ ] `ruff check src tests` is clean
- [ ] New public functions have docstrings explaining their contract
- [ ] New persisted records have a Pydantic model and are JSONL-roundtripped
- [ ] CHANGELOG entry is updated
