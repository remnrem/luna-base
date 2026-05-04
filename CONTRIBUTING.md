# Contributing to Luna

Thank you for contributing. This repository is the core C/C++ implementation of the Luna ecosystem, so changes here can affect the command-line tool, the embeddable library, and downstream interfaces built on top of it.

## Bug reports

Open an issue on [GitHub Issues](https://github.com/remnrem/luna-base/issues) and include:

- operating system and architecture
- compiler and version
- the exact build command you ran
- whether `LGBM=1` or `LGBM=0` was used
- the Luna command, sample-list snippet, or minimal code example that triggers the problem
- the observed output, error text, or crash details

If the issue is data-dependent, reduce it to the smallest reproducible example you can share.

## Feature requests

Open an issue describing:

- the research or engineering use case
- the current limitation
- the proposed command, option, API extension, or behavior change

For command-surface changes, include the expected syntax and output shape if possible.

## Pull requests

1. Fork the repository and create a focused branch from the current default branch.
2. Keep the change set narrow: one logical fix or feature per pull request.
3. Update tests, docs, or help text when the behavior or interface changes.
4. Build the relevant targets before submitting:
   ```bash
   make ARCH=LINUX
   ```
   or
   ```bash
   make ARCH=MAC
   ```
5. Run the built-in test target when practical:
   ```bash
   make test
   ```
6. Open the pull request with a concise summary of the change, the motivation, and any compatibility implications.

## Code and review expectations

- Prefer small, reviewable patches over broad refactors.
- Preserve existing command semantics unless a deliberate interface change is being proposed.
- If you add or rename command options, update the relevant documentation and help definitions.
- If you touch shared engine behavior, consider downstream impact on `lunapi` and Lunascope.

## Style

There is no requirement to match a separate formatter before opening a pull request, but consistency with the surrounding C/C++ style is expected. Favor clear, minimal changes over cosmetic churn.

## Questions

For project questions that are not bug reports, open an issue or contact the maintainers through the main Luna documentation site: https://zzz.nyspi.org/luna/
