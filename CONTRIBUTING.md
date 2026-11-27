# Contributing to `spatial-audit`

Thank you for your interest in contributing to **`spatial-audit`**! This document outlines how to report issues, request features, contribute code or public-cohort adapters, and participate in the community.

## Getting Support & Reporting Issues

- **Bug reports & feature requests**: Open an issue on GitHub with a clear description, your Python and OS versions, and—if possible—a minimal reproducible manifest or command.
- **Questions about methodology or input formats**: Check [docs/INPUTS.md](docs/INPUTS.md), [docs/COHORTS.md](docs/COHORTS.md), and [docs/SOTA.md](docs/SOTA.md) first, or open a GitHub Discussion / Issue labeled `question`.
- **Academic neutrality**: Discussions of published spatial-omics cohorts must remain objective and methodological, focusing on statistical resolution, design identifiability, and confounding sensitivity rather than making accusations about original studies.

## Development Setup

1. Fork and clone the repository in a **Python 3.13** environment.
2. Install the package with development, agent, archive, and browser extras:
   ```sh
   python -m pip install -e ".[dev,agent,archives,browser]"
   python -m playwright install chromium
   ```
3. Run the verification suite before submitting a pull request:
   ```sh
   pytest -q
   ruff check .
   mypy src scripts
   ```
4. For changes affecting reports, schemas, or the dashboard, also run:
   ```sh
   python scripts/verify_release_gate.py --core-only
   ```

## Coding & Provenance Standards

- **Cross-platform paths**: Never hardcode OS-specific drive letters or absolute local paths; use `pathlib.Path` relative to the repository or manifest.
- **Strict JSON**: All JSON outputs must be strictly RFC-8259 compliant (`allow_nan=False`; represent undefined values as `null`).
- **Provenance integrity**: Public archive adapters must record source URLs, byte sizes, and SHA-256 digests. Synthetic benchmarks must be labeled explicitly as `synthetic` and never mixed into `results/public/`.
- **Tests**: New audit dimensions, adapters, or CLI commands should include unit tests in `tests/`.

## Code of Conduct

All contributors and maintainers are expected to foster a welcoming, respectful, and harassment-free environment. Be constructive in reviews, focus criticism on code and statistical assumptions rather than individuals, and respect diverse scientific backgrounds.
