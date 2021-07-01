# -*- makefile -*-
# this makefile contains rules for managing TICA code
# Copyright 2018 MIT
# Licensed under terms of the GPLv3

VERSION=$(shell grep ^__version__ tica/__init__.py | sed -e 's/__version__ *= *//' | sed -e "s/'//g")

TESTDIR=$(TMPDIR)/tica-test-$(USER)
RELDIR=/sw/tica-versions
RELEASE_LABEL=tica-$(VERSION)

all: help

help:
	@echo "options include:"
	@echo "    test       run all of the tests"
	@echo "    test-units run the unit tests"
	@echo "    install    install the software"
	@echo "    release    install to the release space"
	@echo "    bless      make the latest release the blessed release"

info:
	@echo "    VERSION: $(VERSION)"
	@echo "    TESTDIR: $(TESTDIR)"

coverage:
	make test COVER='coverage run -p'
	coverage combine
	coverage html --directory $(TESTDIR)/coverage-report
	@echo "coverage report is at $(TESTDIR)/coverage-report"

testdir:
	mkdir -p $(TESTDIR)

test: test-units

ifdef SUITE
TEST_OPTS=--test-suite $(SUITE)
endif
test-units:
	$(COVER) ./setup.py test $(TEST_OPTS)

# FIXME: prefer to use setup.py and setup.cfg, but that does not work yet
ifndef $PREFIX
PREFIX=/opt/tica
endif
install:
	mkdir $(PREFIX)
	rsync -arv bin $(PREFIX)
	rsync -arv tica $(PREFIX)
#	rsync -arv examples $(PREFIX)
	rsync -arv calibration_* $(PREFIX)
	cp -p README.md LICENSE AUTHORS $(PREFIX)
	find $(PREFIX) -name "*~" -exec rm {} \;
	python -m compileall $(PREFIX)/tica

release:
	make install PREFIX=$(RELDIR)/$(RELEASE_LABEL)

bless:
	rm -f /sw/tica
	ln -s tica-versions/$(RELEASE_LABEL) /sw/tica
