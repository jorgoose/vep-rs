# Agent guide

vep-rs is a Rust port of Ensembl VEP release 116.2. Ensembl VEP's behaviour is the specification;
[CONTRIBUTING.md](CONTRIBUTING.md) is the source of truth for conventions, and this file caches the
parts that fail CI or review when missed.

## Where the truth lives

- **Spec**: Ensembl VEP's Perl source, read from the pinned image
  `ensemblorg/ensembl-vep:release_116.2` (Step 2 of the close-parity-row skill gives the paths).
  Cite what you match as `Module.pm:line`, in the test or comment that depends on it.
- **Backlog**: `tests/parity/matrix.tsv`, one row per Ensembl VEP flag, feature, plugin and cache
  field; [tests/parity/README.md](tests/parity/README.md) defines the columns. To close a row, follow
  [.claude/skills/close-parity-row/SKILL.md](.claude/skills/close-parity-row/SKILL.md).
- **Oracle**: the golden corpora under `tests/golden/` (Ensembl VEP's output, every column of every
  format) and the concordance harness in `scripts/concordance/`.

## Gotchas

- Sign every commit: `git commit -s`. The DCO check fails a pull request with an unsigned commit.
- A new `.rs` file opens with the two-line license header, a blank line, then its `//!` module doc
  (`scripts/check_module_docs.sh` checks the blank line).
- A flag added or changed in `crates/vep-cli/src/args.rs` gets its matrix row in the same change, and
  the generated CLI reference is rewritten and committed:
  `VEP_UPDATE_CLI_REFERENCE=1 cargo test -p vep-cli --test cli_reference`.
- A comment is one or two present-tense sentences stating a fact the code cannot, usually the Ensembl
  VEP behaviour matched, cited. It names nothing outside the code.
- Diagnostics go through `tracing`; a refusal ends the run with one `ERROR: [<component>] <message>`
  line. CONTRIBUTING.md's Code Conventions lists the few direct writes.
- `CHANGELOG.md` entries follow the grammar in CONTRIBUTING.md's Changelog entries. The trailing
  `(#N)` is added in a commit after the pull request exists.
- A change that can alter annotation output runs `cargo test -p vep-cli --test golden --test
  format_parity`. A consequence or plugin change also runs the concordance smoke (`just
  concordance-smoke <perl_cache> <json_cache> <fasta>`), which needs Docker and an Ensembl VEP 116
  cache; on a host without them, the pull request's Verification section says the smoke was not run
  and why.

## Before pushing

Run, and quote the result in the pull request's Verification section:

```bash
cargo fmt --all -- --check && cargo clippy --workspace -- -D warnings
cargo test --workspace
python3 scripts/parity/check_parity_matrix.py && scripts/check_module_docs.sh
```

Fill in [.github/pull_request_template.md](.github/pull_request_template.md) for the description.

## Scope

- One matrix row or one issue per branch, cut from `origin/main`.
- A row's `target` is the maintainers' call. When the Perl source or the code shows the target is
  wrong, stop and report it as a **target decision**, with the evidence, and keep the row as it is.
- Parallel branches collide on `tests/parity/matrix.tsv`, `CHANGELOG.md`,
  `crates/vep-cli/src/args.rs`, `docs/cli-reference.md` and `docs/cli/`. Rebase onto `origin/main`
  immediately before pushing and regenerate the CLI reference after the rebase.
- Releases stay with the maintainers: the version in `Cargo.toml`, tags, crates.io publishing, and the
  `release.yml`, `docker.yml` and `pages.yml` workflows.
