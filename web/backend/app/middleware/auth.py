"""Auth middleware — pass-through sentinel for future implementation."""

from fastapi import Request
from starlette.middleware.base import BaseHTTPMiddleware


class AuthMiddleware(BaseHTTPMiddleware):
    """Placeholder auth middleware. Currently passes all requests through."""

    async def dispatch(self, request: Request, call_next):
        # Future: validate API keys, JWT tokens, etc.
        response = await call_next(request)
        return response
