# Development guide

## Setup

Requires Python 3.10 or newer (3.11 recommended).

    git clone https://github.com/gmrandazzo/PAutoDock.git
    cd PAutoDock
    python3.11 -m venv .venv
    source .venv/bin/activate
    pip install -e . pytest pytest-cov flake8 tox pre-commit mypy

or with Poetry:

    poetry install

## Running the tests

    pytest tests

or through tox, which runs the test suite and the type check:

    tox -e py3.11,mypy

The test suite mocks the external binaries (obabel, vina, autodock4),
so no docking software is needed to develop.

## Type checking

The codebase is fully annotated and checked with mypy in strict mode:

    mypy src

The configuration lives in `pyproject.toml` (`[tool.mypy]`).

## Linting and formatting

Run the whole pre-commit suite (autopep8, black, isort, flake8, yesqa,
mypy and the standard file checks):

    pre-commit run --all-files

## Project layout

- `src/pautodock/__main__.py` -- command line interface
- `src/pautodock/adparallel.py` -- the screening pipeline (preparation,
  parallel execution, results collection)
- `src/pautodock/molop.py` -- molecule and receptor conversions and
  coordinate helpers (Open Babel based)
- `src/pautodock/multimol2op.py` -- multi-mol2 splitting
- `src/pautodock/fileutils.py` -- executable lookup and archive helpers
- `src/pautodock/mgltoolsinstall.py` -- optional MGLTools installer
- `src/pautodock/__recover_output__.py` -- results recovery tool
- `src/pautodock/__autogridmap2dx__.py` -- AutoGrid to OpenDX converter

## Branching model

Development happens on the `develop` branch; `main` always holds the
latest release.

1. Work on `develop` and open a pull request to `main`. Every pull
   request and every push runs the CI: the test suite, the mypy type
   check and the pre-commit suite (`.github/workflows/pytest.yml` and
   `.github/workflows/pre-commit.yml`).
2. Merge the pull request.

## Release cycle

1. Bump the version in both `pyproject.toml` and
   `src/pautodock/__init__.py`.
2. Work on `develop` and merge it into `main` through a pull request.
3. Tag the merge commit on `main` and push the tag:

       git fetch origin
       git tag v<X.Y.Z> origin/main
       git push origin v<X.Y.Z>

Pushing the tag runs the release workflow
(`.github/workflows/release.yml`): the tests run first, then the
package is built with Poetry, a GitHub Release is created with the
source distribution and the wheel attached, and the package is
published to PyPI with twine.

The PyPI upload uses the `PYPI_API_TOKEN` repository secret (a PyPI API
token; configure under Settings -> Secrets and variables -> Actions).
The GitHub Release uses the automatic `GITHUB_TOKEN`.
