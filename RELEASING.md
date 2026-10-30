# Releasing Parachron on PyPI

Parachron is published manually, without CI, and **only as a source distribution (sdist)**: pip compiles it on
install, which needs a C compiler. No binary wheels are published (decision of Oct 30, 2026): the speed of the
WHFast512 kernels depends on the compiler (GCC 13+ recommended, clang 9–13 % slower), and the README tells users how
to select it; a source build also cannot carry a compiler regression of a wheel build into every installation.

## One-time setup

```bash
python3 -m pip install --upgrade build twine
```

Create an API token on https://pypi.org (Account settings → API tokens). For a first upload it must be scoped to the
whole account, because the project `parachron` does not exist on PyPI yet. After the first upload, replace it by a
token scoped to `parachron`. Optionally do the same on https://test.pypi.org for a trial run.

## Release

1. Set the version in `pyproject.toml` (`version = "..."`); this is the package version on PyPI. The C library keeps
   reporting the REBOUND base version (`parachron.__version__`).
2. Build from a clean checkout of the release commit:
   ```bash
   git clone https://github.com/hippke/parachron && cd parachron
   python3 -m build --sdist
   python3 -m twine check dist/*
   ```
3. Optional trial on TestPyPI:
   ```bash
   python3 -m twine upload --repository testpypi dist/*
   python3 -m venv /tmp/pt && /tmp/pt/bin/pip install --index-url https://test.pypi.org/simple/ parachron
   /tmp/pt/bin/python -c "import parachron; print(parachron.__version__, parachron.MEP)"
   ```
4. Upload (username `__token__`, password = the API token):
   ```bash
   python3 -m twine upload dist/*
   ```
5. Tag the release: `git tag v<version> && git push origin v<version>`.

## Checks before a release

```bash
pip install -e . && python3 -m pytest -q parachron/tests -k "not horizons and not plotting and not test_method"
```

The source distribution must contain `LICENSE`, `CHANGES_FROM_REBOUND.md` (GPL-3.0 section 5(a)) and
`src/integrator_whfast512.s`. Check with `tar tzf dist/*.tar.gz`.
