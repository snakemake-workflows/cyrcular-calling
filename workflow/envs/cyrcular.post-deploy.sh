# this is needed to have the `ruget` binary installed by `plotly_kaleido`
# available in the PATH, it is removed from PATH after installation of cyrcular
PATH="$HOME/.cargo/bin:$PATH"
# TODO: merge branch `annotate` in cyrcular and create bioconda recipe
cargo install --git https://github.com/tedil/cyrcular.git --branch fix--make-table-command-work-with-any-set-of-event-definitions --rev e6f18dd287ae9233dbcf6fd7d0b05039512816f8 --locked --root "${CONDA_PREFIX}"
PATH=${PATH#[^:]*:}