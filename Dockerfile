FROM python:3.12-slim

WORKDIR /app

# Install uv for fast package installation
COPY --from=ghcr.io/astral-sh/uv:0.4.30 /uv /usr/local/bin/uv

# Copy package files
COPY pyproject.toml uv.lock README.rst LICENSE MANIFEST.in ./
COPY breizorro/ breizorro/

# Install git so uv can pull the scabha repository
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    && rm -rf /var/lib/apt/lists/*
# Install package with all dependencies using the locked uv environment
# Run a locked sync, then install optional extras via pip from the local project.
RUN uv sync --frozen --extra all
# Ensure the venv-installed CLI is available
ENV PATH="/app/.venv/bin:${PATH}"

# Make CLI available
CMD ["breizorro", "--help"]
