#!/bin/bash
# Chromatiq Build Script
# Builds WASM binaries and prepares web distribution

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "🔧 Chromatiq Build Script"
echo "========================"

# Check for wasm-pack
if ! command -v wasm-pack &> /dev/null; then
    echo "📦 Installing wasm-pack..."
    cargo install wasm-pack
fi

# Build WASM module
echo "🏗️  Building Chromatiq WASM module..."
cd crates/chromatiq-wasm
wasm-pack build --target web --out-dir www/pkg

echo "✅ Build complete!"
echo ""
echo "📂 Output: crates/chromatiq-wasm/www/pkg/"
echo ""
echo "🌐 To run the web interface:"
echo "   cd crates/chromatiq-wasm/www && npx serve ."
