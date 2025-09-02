# Use official Python base image
FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Copy all project files into the container
COPY . /app

# Install dependencies
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r requirements.txt

# Optional: Set environment variables for audit/debug
ENV PYTHONUNBUFFERED=1 \
    CSI_ENV=production

# Expose FastAPI port
EXPOSE 8000

# Default command: Launch FastAPI app
ENV PYTHONPATH=/app
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
