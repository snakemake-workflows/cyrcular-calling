echo "PATH: $PATH"
echo "HOME: $HOME"
echo "GITHUB_HOME: $GITHUB_HOME"
echo "GITHUB_PATH: $GITHUB_PATH"
PATH="$HOME/.cargo/bin:$PATH"
cargo install --git https://github.com/tedil/cyrcular.git --branch annotate --rev 44ba5274b9f2aee7d193b774343acaccf0edec53 --locked
