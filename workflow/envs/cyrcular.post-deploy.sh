# this is needed to have the `ruget` binary installed by `plotly_kaleido`
# available in the PATH, it is removed from PATH after installation of cyrcular
PATH="$HOME/.cargo/bin:$PATH"
# TODO: merge branch `annotate` in cyrcular and create bioconda recipe
cargo install --git https://github.com/tedil/cyrcular.git --branch perf/reduce-memory-footprint --rev 46b379797e73099921282c0e6d5ba8ac7d8f0497 --locked --root "${CONDA_PREFIX}"
PATH=${PATH#[^:]*:}