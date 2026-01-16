use clap::Parser;
use tracing_subscriber::EnvFilter;

mod catalog;
mod cli;
mod core;
mod matching;
mod parsing;
mod utils;
mod web;

fn main() -> anyhow::Result<()> {
    let cli = cli::Cli::parse();

    // Initialize logging based on verbosity flag
    let filter = if cli.verbose {
        EnvFilter::new("ref_solver=debug,info")
    } else {
        EnvFilter::new("ref_solver=warn")
    };

    tracing_subscriber::fmt()
        .with_env_filter(filter)
        .with_target(false)
        .without_time()
        .init();

    match cli.command {
        cli::Commands::Identify(args) => {
            cli::identify::run(args, cli.format, cli.verbose)?;
        }
        cli::Commands::Compare(args) => {
            cli::compare::run(args, cli.format, cli.verbose)?;
        }
        cli::Commands::Score(args) => {
            cli::score::run(args, cli.format, cli.verbose)?;
        }
        cli::Commands::Catalog(args) => {
            cli::catalog::run(args, cli.format, cli.verbose)?;
        }
        cli::Commands::Serve(args) => {
            web::server::run(args)?;
        }
    }

    Ok(())
}
