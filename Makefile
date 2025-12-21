# Build the Docker image
build:
	docker build -t dps_pipeline .

# Run the pipeline using Docker
run:
	docker run dps_pipeline

# Build + run in one command
all: build run