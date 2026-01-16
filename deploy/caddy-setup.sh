#!/bin/bash
# Caddy setup script for ref-solver on Lightsail
# Run this on your Lightsail instance:
#   curl -sSL https://raw.githubusercontent.com/fulcrumgenomics/ref-solver/main/deploy/caddy-setup.sh | bash

set -e

echo "=== Setting up ref-solver with Caddy ==="

# Install Caddy
echo "Installing Caddy..."
sudo apt-get update
sudo apt-get install -y debian-keyring debian-archive-keyring apt-transport-https curl
curl -1sLf 'https://dl.cloudsmith.io/public/caddy/stable/gpg.key' | sudo gpg --dearmor -o /usr/share/keyrings/caddy-stable-archive-keyring.gpg
curl -1sLf 'https://dl.cloudsmith.io/public/caddy/stable/debian.deb.txt' | sudo tee /etc/apt/sources.list.d/caddy-stable.list
sudo apt-get update
sudo apt-get install -y caddy

# Stop nginx if running (we're replacing it)
sudo systemctl stop nginx 2>/dev/null || true
sudo systemctl disable nginx 2>/dev/null || true

# Create ref-solver service (bind to localhost since Caddy will proxy)
echo "Configuring ref-solver service..."
sudo tee /etc/systemd/system/ref-solver.service > /dev/null <<'EOF'
[Unit]
Description=ref-solver web server
After=network.target

[Service]
Type=simple
User=www-data
ExecStart=/usr/local/bin/ref-solver serve --address 127.0.0.1 --port 8080
Restart=always
RestartSec=5
Environment=RUST_LOG=info

[Install]
WantedBy=multi-user.target
EOF

# Configure Caddy
echo "Configuring Caddy..."
sudo tee /etc/caddy/Caddyfile > /dev/null <<'EOF'
# Primary domain
refsolver.bio, www.refsolver.bio {
    reverse_proxy localhost:8080
}

# Redirect .com to .bio
refsolver.com, www.refsolver.com {
    redir https://refsolver.bio{uri} permanent
}
EOF

# Reload systemd and start services
echo "Starting services..."
sudo systemctl daemon-reload
sudo systemctl enable ref-solver
sudo systemctl restart ref-solver
sudo systemctl enable caddy
sudo systemctl restart caddy

echo ""
echo "=== Setup complete ==="
echo ""
echo "Caddy will automatically obtain SSL certificates from Let's Encrypt."
echo ""
echo "Your sites:"
echo "  https://refsolver.bio     (primary)"
echo "  https://www.refsolver.bio (primary)"
echo "  https://refsolver.com     -> redirects to refsolver.bio"
echo "  https://www.refsolver.com -> redirects to refsolver.bio"
echo ""
echo "Check status with:"
echo "  sudo systemctl status caddy"
echo "  sudo systemctl status ref-solver"
echo ""
