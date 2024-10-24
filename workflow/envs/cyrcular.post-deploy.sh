# this is needed to have the `ruget` binary installed by `plotly_kaleido`
# available in the PATH, it is removed from PATH after installation of cyrcular
PATH="$HOME/.cargo/bin:$PATH"
# TODO: merge branch `annotate` in cyrcular and create bioconda recipe
cargo install --git https://github.com/tedil/cyrcular.git --branch fix--make-table-command-work-with-any-set-of-event-definitions --rev ac7a6d8d3bdb15a57d621e1d7f588f6ee875ae91 --locked --root "${CONDA_PREFIX}"
PATH=${PATH#[^:]*:}