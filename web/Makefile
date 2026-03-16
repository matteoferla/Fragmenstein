.PHONY: dev backend frontend install clean

# Start both backend and frontend for development
dev: backend frontend

backend:
	cd backend && uvicorn app.main:app --reload --port 8000

frontend:
	cd frontend && npm run dev

install:
	cd backend && pip install -r requirements.txt
	cd frontend && npm install

clean:
	rm -rf backend/data
	rm -rf frontend/.next frontend/node_modules
