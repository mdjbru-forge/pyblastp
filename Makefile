### * Variables

PYTHON_MODULE=pyblastp
PYTHON_MODULE_EGG=$(PYTHON_MODULE).egg-info

COVERED_PACKAGES=$(PYTHON_MODULE)
SPHINX_DOC_FOLDER=docs/

### * Help

help:
	@echo "Makefile for the $(PYTHON_MODULE) Python module                "
	@echo "                                                               "
	@echo "Type \"make <target>\" where <target> is one of the following: "
	@echo "                                                               "
	@echo "  test        Run the tests with coverage output               "
	@echo "  doc         Run sphinx to generate the module documentation  "
	@echo "  clean       Clean everything (coverage, docs, pyc files)     "
	@echo "                                                               "
	@echo "  install     Install the module and command-line tools        "
	@echo "  uninstall   Uninstall the module                             "

### * Main targets

### ** test
tests: test
test:
	nosetests tests/ --with-coverage --cover-package=$(COVERED_PACKAGES) --cover-html \
          --with-html --html-file=tests/nosetests.html
	@echo -e "\nThe coverage results are accessible from cover/index.html"
	@echo "The html version of the test results are accessible from tests/nosetests.html"

### ** doc
doc: moduleGraph
	sphinx-apidoc -o $(SPHINX_DOC_FOLDER) $(PYTHON_MODULE)/
	cd $(SPHINX_DOC_FOLDER); make html
	@echo -e "\nThe documentation is accessible from $(SPHINX_DOC_FOLDER)_build/html/index.html"

### ** moduleGraph
moduleGraph :
	rm -f $(SPHINX_DOC_FOLDER)/*.sphinx.dot
	makeSphinxGraphs $(PYTHON_MODULE)/*.py -o $(SPHINX_DOC_FOLDER)

### ** clean
clean:
	rm -f .coverage
	rm -fr cover
	rm -f tests/nosetests.html
	rm -fr docs/_build
	# http://superuser.com/questions/112078/delete-matching-files-in-all-subdirectories
	find . -name \*.pyc -type f -delete

### ** install
install:
	rm -fr $(PYTHON_MODULE_EGG)
	pip install --user -e .

### ** uninstall
uninstall:
	rm ~/.local/bin/pyblastp
	pip uninstall -y $(PYTHON_MODULE)
