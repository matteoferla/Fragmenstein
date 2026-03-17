.PHONY: start dev backend frontend install clean

# Single command to start everything
start:
	./start.sh

# Development (same as start)
dev: start

# Start backend only
backend:
	cd backend && python -m uvicorn app.main:app --reload --port 8000

# Start frontend only
frontend:
	cd frontend && npm run dev

# Install all dependencies
install:
	@echo "Installing backend dependencies..."
	@command -v uv >/dev/null && uv pip install --system -r backend/requirements.txt || pip install -r backend/requirements.txt
	@echo "Installing frontend dependencies..."
	cd frontend && npm install
	@echo "Done."

# Remove runtime data (keeps code)
clean:
	rm -rf backend/data backend/output
	rm -rf frontend/.next
