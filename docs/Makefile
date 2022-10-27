# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = raven
SOURCEDIR     = source
BUILDDIR      = build

GHBUILDDIR    = ../../raven-docs


# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O) -W

.PHONY: gh-pages
gh-pages: html
	cp -r $(BUILDDIR)/html/* $(GHBUILDDIR)/html/
	cd $(GHBUILDDIR)/html; git add . ; git commit -m "rebuilt docs"; git push origin gh-pages
	@echo
	@echo "Committed to gh-pages."
