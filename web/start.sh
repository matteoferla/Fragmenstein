#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────
# Fragmenstein Web UI — Single-command launcher
# Usage:  cd web && ./start.sh
# ──────────────────────────────────────────────────────────────

set -e

DIR="$(cd "$(dirname "$0")" && pwd)"
VERSION=$(cat "$DIR/VERSION" 2>/dev/null || echo "0.0.0")

echo ""
echo "  ╔═══════════════════════════════════════════╗"
echo "  ║        Fragmenstein Web UI v${VERSION}         ║"
echo "  ║   Fragment-Based Drug Design Pipeline     ║"
echo "  ╚═══════════════════════════════════════════╝"
echo ""

# ── Check dependencies ────────────────────────────────────────

check_cmd() {
  if ! command -v "$1" &>/dev/null; then
    echo "  ✗ $1 not found. $2"
    return 1
  fi
  return 0
}

OK=true
check_cmd python "Install Python 3.10+" || OK=false
check_cmd node "Install Node.js 18+" || OK=false
check_cmd npm "Comes with Node.js" || OK=false

if [ "$OK" = false ]; then
  echo ""
  echo "  Please install missing dependencies and try again."
  exit 1
fi

# ── Install if needed ─────────────────────────────────────────

if [ ! -d "$DIR/frontend/node_modules" ]; then
  echo "  Installing frontend dependencies..."
  (cd "$DIR/frontend" && npm install --silent)
  echo "  ✓ Frontend dependencies installed"
fi

# Check Python deps
python -c "import fastapi, uvicorn, rdkit" 2>/dev/null || {
  echo "  Installing backend dependencies..."
  if command -v uv &>/dev/null; then
    uv pip install --system -r "$DIR/backend/requirements.txt" --quiet
  else
    pip install -r "$DIR/backend/requirements.txt" --quiet
  fi
  echo "  ✓ Backend dependencies installed"
}

# ── Launch ────────────────────────────────────────────────────

BACKEND_PORT=${BACKEND_PORT:-8000}
FRONTEND_PORT=${FRONTEND_PORT:-3000}

echo "  Starting backend on port ${BACKEND_PORT}..."
(cd "$DIR/backend" && python -m uvicorn app.main:app --host 0.0.0.0 --port "$BACKEND_PORT" --reload) &
BACKEND_PID=$!

echo "  Starting frontend on port ${FRONTEND_PORT}..."
(cd "$DIR/frontend" && NEXT_PUBLIC_API_URL="http://localhost:${BACKEND_PORT}" npm run dev -- -p "$FRONTEND_PORT") &
FRONTEND_PID=$!

# ── Wait for startup ──────────────────────────────────────────

sleep 3
echo ""
echo "  ┌───────────────────────────────────────────┐"
echo "  │                                           │"
echo "  │   Open in your browser:                   │"
echo "  │   http://localhost:${FRONTEND_PORT}                   │"
echo "  │                                           │"
echo "  │   Press Ctrl+C to stop                    │"
echo "  │                                           │"
echo "  └───────────────────────────────────────────┘"
echo ""

# ── Handle shutdown ───────────────────────────────────────────

cleanup() {
  echo ""
  echo "  Shutting down..."
  kill "$BACKEND_PID" 2>/dev/null
  kill "$FRONTEND_PID" 2>/dev/null
  wait "$BACKEND_PID" 2>/dev/null
  wait "$FRONTEND_PID" 2>/dev/null
  echo "  Done."
}

trap cleanup EXIT INT TERM

# Wait for either process to exit
wait
