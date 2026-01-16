#!/bin/bash
# Run this on your Lightsail instance to set up the server
# Usage: curl -sSL <raw-github-url> | bash

set -e

echo "=== Setting up ref-solver server ==="

# Update system
sudo apt-get update
sudo apt-get install -y nginx certbot python3-certbot-nginx

# Create service file
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

# Enable service (will start after first deploy)
sudo systemctl daemon-reload
sudo systemctl enable ref-solver

# Configure nginx as reverse proxy
sudo tee /etc/nginx/sites-available/ref-solver > /dev/null <<'EOF'
server {
    listen 80;
    server_name _;

    location / {
        proxy_pass http://127.0.0.1:8080;
        proxy_http_version 1.1;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;

        # For file uploads
        client_max_body_size 20M;
    }
}
EOF

sudo ln -sf /etc/nginx/sites-available/ref-solver /etc/nginx/sites-enabled/
sudo rm -f /etc/nginx/sites-enabled/default
sudo nginx -t && sudo systemctl reload nginx

echo ""
echo "=== Setup complete ==="
echo ""
echo "Next steps:"
echo "1. Deploy the binary via GitHub Actions (push to main)"
echo "2. For HTTPS with a domain, run:"
echo "   sudo certbot --nginx -d yourdomain.com"
echo ""
