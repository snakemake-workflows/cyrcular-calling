//! This is a regular crate doc comment, but it also contains a partial
//! Cargo manifest.  Note the use of a *fenced* code block, and the
//! `cargo` "language".
//!
//! ```cargo
//! [dependencies]
//! csv = "1.1"
//! rust-bio = "*"
//! ```


use anyhow::Result;

fn main() -> Result<()> {
    let _stdout_redirect = snakemake.redirect_stdout(snakemake.log.stdout)?;
    println!("This will be written to path/to/stdout.log");

    // redirect stderr to "path/to/stderr.log"
    let _stderr_redirect = snakemake.redirect_stderr(snakemake.log.stderr)?;
    Ok(())
}