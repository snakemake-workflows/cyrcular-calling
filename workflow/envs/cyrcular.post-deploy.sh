# this is needed to have the `ruget` installed by `plotly_kaleido`
# available in the PATH, it is removed from PATH after installed
PATH="$HOME/.cargo/bin:$PATH"
cargo install --git https://github.com/tedil/cyrcular.git --branch annotate --rev 44ba5274b9f2aee7d193b774343acaccf0edec53 --locked --root "${CONDA_PREFIX}" --no-track
PATH=${PATH#[^:]*:}
