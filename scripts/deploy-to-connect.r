#!/usr/bin/env Rscript

# Deploy to Posit Connect

cat("Starting deployment script...\n")

# Load environment variables from .env if present
if (file.exists(".env")) {
  suppressPackageStartupMessages(library(dotenv))
  load_dot_env()
  cat("Loaded .env file\n")
}

# Load rsconnect
cat("Loading rsconnect package...\n")
suppressPackageStartupMessages(library(rsconnect))

# Get credentials
server_url <- Sys.getenv("CONNECT_SERVER")
api_key <- Sys.getenv("CONNECT_API_KEY")
account <- Sys.getenv("CONNECT_ACCOUNT")

cat("Server URL:", if(server_url == "") "MISSING" else "SET", "\n")
cat("API Key:", if(api_key == "") "MISSING" else "SET", "\n")
cat("Account:", if(account == "") "MISSING" else "SET", "\n")

# Validate credentials
if (any(c(server_url, api_key, account) == "")) {
  stop("Missing required environment variables: CONNECT_SERVER, CONNECT_API_KEY, or CONNECT_ACCOUNT")
}

# Check output directory
if (!dir.exists("_output")) {
  stop("Output directory '_output' not found. Run quarto render first.")
}
cat("Output directory exists\n")

# Extract server name and ensure rsconnect directories exist
server_name <- gsub("https?://", "", server_url)
cat("Server name:", server_name, "\n")

cat("Deploying to", server_name, "...\n")

# Create rsconnect config directories
rsconnect_dir <- file.path(Sys.getenv("HOME"), ".config", "R", "rsconnect")
servers_dir <- file.path(rsconnect_dir, "servers")
accounts_dir <- file.path(rsconnect_dir, "accounts")

cat("Creating rsconnect directories...\n")
dir.create(rsconnect_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(servers_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(accounts_dir, recursive = TRUE, showWarnings = FALSE)

# Add server and connect
cat("Adding server...\n")
tryCatch({
  addServer(url = server_url, name = server_name)
  cat("Server added successfully\n")
}, error = function(e) {
  cat("Error adding server:", as.character(e), "\n")
  stop(e)
})

cat("Connecting API user...\n")
tryCatch({
  connectApiUser(account = account, server = server_name, apiKey = api_key)
  cat("API user connected successfully\n")
}, error = function(e) {
  cat("Error connecting API user:", as.character(e), "\n")
  stop(e)
})

# Set a specific app ID to ensure we update the same app
app_id <- "721"  # Use the ID from your first successful deployment
cat("Using app ID:", app_id, "\n")

# Deploy as static site with specific app ID
cat("Starting deployment...\n")
tryCatch({
  deployApp(
    appDir = "_output",
    appId = app_id,
    account = account,
    server = server_name,
    forceUpdate = TRUE,
    contentCategory = "site"
  )
  cat("Deployment completed successfully!\n")
}, error = function(e) {
  cat("Error during deployment:", as.character(e), "\n")
  stop(e)
})

cat("Script finished!\n")