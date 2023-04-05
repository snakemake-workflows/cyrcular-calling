echo "CARGO_HOME: $CARGO_HOME"
echo "$CARGO_HOME/bin" >> PATH
cargo install --git https://github.com/tedil/cyrcular.git --branch annotate --rev 44ba5274b9f2aee7d193b774343acaccf0edec53 --locked
