# Luna

![lunabase CI](https://github.com/remnrem/luna-base/workflows/lunabase%20CI/badge.svg) [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**Luna is the core C/C++ analysis engine of the Luna sleep research ecosystem.** It provides command-line and embeddable tooling for large-scale analysis of EDF/EDF+ polysomnography data, annotation handling, signal processing, sleep staging, spectral analysis, cohort summaries, and derived-model workflows.  Please see the [primary documentation](https://zzz.nyspi.org/luna/) for more details on usage and installation.

![Lunascope interface built on Luna](docs/assets/lunascope-interface.webp)

## Technical report

For a higher-level overview of how Luna fits together with `lunapi` and Lunascope, see the technical report, [*The Luna ecosystem: integrated tools for sleep signal visualization and analysis*](https://doi.org/10.5281/zenodo.20025864). It describes the broader Luna ecosystem and the shared workflow spanning scripted and graphical analysis.

## What Luna provides

Luna is designed for reproducible sleep research workflows that need to scale from a single recording to large cohorts while preserving a common analytical core across command-line, C++ embedding, Python, and interactive desktop review.

- Native C/C++ tooling for EDF/EDF+, annotations, masking, restructuring, and batch execution
- A broad command set for spectral analysis, spindle and slow-wave detection, connectivity, artifact handling, actigraphy, staging, and predictive models
- A reusable static/shared library (`libluna.a`, `libluna.dylib`, `libluna.so`) for embedding Luna in other applications
- A command-line executable (`luna`) for scripted and cohort-scale processing
- Shared analytical infrastructure used by `lunapi` and Lunascope

## Ecosystem

Luna is one part of a larger stack. The same analytical engine is exposed through multiple interfaces depending on the workflow:

| Project | Role | Repository | Web documentation |
|---|---|---|---|
| `luna-base` | Core C/C++ engine and CLI | https://github.com/remnrem/luna-base | https://zzz.nyspi.org/luna/ |
| `lunapi` | Python interface for notebooks and scripted workflows | https://github.com/remnrem/luna-api | https://zzz.nyspi.org/luna/python/ |
| `Lunascope` | Interactive desktop application for visualization, annotation, and review | https://github.com/Lorcan7274/lunascope | https://zzz.nyspi.org/luna/lunascope/ |

## Quick usage

Luna processes either a single EDF or a [sample list](https://zzz.nyspi.org/luna/luna/args/#sample-lists) plus one or more commands:

```bash
./luna sample.lst -o out.db -s 'EPOCH & PSD sig=C3 spectrum'
```

### Useful built-in help:

```bash
./luna -h
./luna -h PSD
./luna -hs spindle
```

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for bug reports, feature requests, and pull request guidance specific to this C++ codebase.
