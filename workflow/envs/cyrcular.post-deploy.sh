# this is needed to have the `ruget` installed by `plotly_kaleido`
# available in the PATH, it is removed from PATH after installed
PATH="$HOME/.cargo/bin:$PATH"
PATH=${PATH#[^:]*:}
