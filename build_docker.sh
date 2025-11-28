#!/bin/bash

# Extract version from config.py
VERSION=$(grep "__version__" src/config.py | cut -d "'" -f 2)

echo "Building SQANTI3 Docker image version ${VERSION}..."

# Build Docker image with version as build argument
docker build \
    --no-cache \
    --build-arg SQANTI3_VERSION="${VERSION}" \
    -t anaconesalab/sqanti3:v${VERSION} \
    -t anaconesalab/sqanti3:latest \
    .
# Only if build is successful, print success message
if [ $? -ne 0 ]; then
    echo "Docker build failed!"
    exit 1
fi
echo "Docker image anaconesalab/sqanti3:${VERSION} built successfully!"