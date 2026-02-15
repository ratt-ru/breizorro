FROM python:3.12-slim

WORKDIR /app

# Install uv for fast package installation
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

# Copy package files
COPY pyproject.toml README.rst LICENSE MANIFEST.in ./
COPY breizorro/ breizorro/

# Install package with all dependencies using uv (much faster than pip)
# Install with [all] extras to include catalog and gui functionality
RUN uv pip install --system --no-cache ".[all]"

# Make CLI available
CMD ["breizorro", "--help"]
