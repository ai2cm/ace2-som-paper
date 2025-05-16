ENVIRONMENT_NAME=ace2-som-paper-experiments
BASE_ENVIRONMENT=base-environment.yaml
UNCOMPILED_REQUIREMENTS=requirements.in
COMPILED_REQUIREMENTS=requirements.txt


create_base_environment:
	conda env create --name $(ENVIRONMENT_NAME) --file $(BASE_ENVIRONMENT)


lock_pip:
	conda run --no-capture-output -n $(ENVIRONMENT_NAME) \
		pip-compile --strip-extras $(UNCOMPILED_REQUIREMENTS) > $(COMPILED_REQUIREMENTS)


install_pip_dependencies:
	conda run --no-capture-output -n $(ENVIRONMENT_NAME) pip install -r $(COMPILED_REQUIREMENTS)


create_environment: create_base_environment lock_pip install_pip_dependencies
