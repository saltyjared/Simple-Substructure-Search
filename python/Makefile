# Makefile for creating a Python environment and running tests for the substructure searches

# Create a Python environment in conda with all necessary packages
environment: environment.yaml
	conda env create -f environment.yaml

# Run all associated test cases
test:
	pytest -v