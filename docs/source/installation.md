# Installation

## Standard Installation

This module is distributed using [`PyPI`](https://pypi.org/project/idea-bio/)
under the name `idea-bio`.

We can install it using python's package manager `pip`:

```bash
pip install idea-bio
```

## Additional Dependency: `cargo`

If you were able to install the module using the `pip` command above
then feel free to move on, otherwise if you received an error asking
for `cargo` in your `$PATH` then this is for you!

This module relies on [`ggetrs`](https://github.com/noamteyssier/ggetrs)
for mediating the gene set enrichment analysis, which is built in rust.

As a result we require `cargo`, the rust package manager to include it as
a dependency.
We can install `cargo` easily with the following line:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

We then can either reload our terminal (close and reopen) or source
the environment to put `cargo` into our `$PATH`:

```bash
source ~/.cargo/env
```

Finally we can install via `pip`:

```bash
pip install idea-bio
```
