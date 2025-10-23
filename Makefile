# Makefile for PALM-meteo documentation and release

# Copyright 2018-2025 Institute of Computer Science of the Czech Academy of
# Sciences, Prague, Czech Republic. Authors: Pavel Krc
#
# This file is part of PALM-METEO.
#
# PALM-METEO is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# PALM-METEO is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# PALM-METEO. If not, see <https://www.gnu.org/licenses/>.

DOXYGEN = doxygen
BROWSER = firefox
RUNNER  = python3 -m palmmeteo
TWINE   = python3 -m twine

SOURCES_MAIN = palmmeteo/*.py palmmeteo/config_init.yaml
SOURCES_PLG  = palmmeteo_stdplugins/*.py palmmeteo_stdplugins/config_init.yaml
SOURCES      = $(SOURCES_MAIN) $(SOURCES_PLG)
DOXYCONFIG   = docs/doxygen.config
DOCS_ROOT    = docs/html
DOCS_INDEX   = $(DOCS_ROOT)/index.html
DOCS_PAGES   = docs/pages/*
CMDLINE_DOC  = docs/examples/commandline.txt
EXAMPLES     = $(CMDLINE_DOC) docs/examples/template.yaml

# === project building ===

build: dist/.updated

dist/.updated: $(SOURCES) $(DOCS_INDEX)
	python3 -m build -n && touch $@

upload: dist/.uploaded

dist/.uploaded: dist/*
	twine upload $? && touch $@

# === docs ===

docs_build: $(DOCS_INDEX)

$(DOCS_INDEX): $(DOXYCONFIG) $(SOURCES) README.md $(DOCS_PAGES) $(EXAMPLES)
	set -e
	rm -rf docs/html/*
	$(DOXYGEN) $(DOXYCONFIG)
	touch $(DOCS_ROOT)/.nojekyll
	cd docs/html && [ gh-pages = "`git branch --show-current`" ] && git add -A && git commit -m "Update github pages"

$(CMDLINE_DOC): $(SOURCES_MAIN)
	$(RUNNER) -h > $@

show: $(DOCS_INDEX)
	$(BROWSER) $(DOCS_INDEX)

docs_publish: $(DOCS_INDEX)
	git push github gh-pages

all: build upload

clean:
	rm -rf docs/html docs/man docs/latex $(CMDLINE_DOC) dist/*

.PHONY: build upload docs_build show docs_publish all clean

# Notes for using the subdirectory docs/html as a separate worktree for
# the orphan branch gh-pages:
# After fetching gh-pages from github, just run
# 	git worktree add docs/html gh-pages
# Occassionally clean up history (from docs/html):
# 	git checkout github/gh-pages
#	git branch -d gh-pages
#	git checkout --orphan gh-pages
#	+ add all, commit, force push
#
# Previously tried using docs/html as a submodule (also works, but needs more
# maintenance + storage) and as a pushed subtree (even worse).
