# Epilect developer Makefile
PYTHON ?= python

.PHONY: install test demo clean help

help:
	@echo "make install   - pip install -e . (editable, with dev extras)"
	@echo "make test      - run the pytest unit suite"
	@echo "make demo      - classify the bundled SDSE mock and regenerate demo/demo_panel.png"
	@echo "make clean     - remove build/test/demo artifacts"

install:
	$(PYTHON) -m pip install -e ".[dev]"

test:
	$(PYTHON) -m pytest -q

demo:
	$(PYTHON) demo/run_demo.py

clean:
	rm -rf build dist *.egg-info .pytest_cache
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
	rm -f demo/demo_panel.png demo/demo_classify.tsv
