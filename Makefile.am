ACLOCAL_AMFLAGS = -I config/m4 -I config/m4/acarchive

SUBDIRS=

SUBDIRS += scripts
SUBDIRS += src
SUBDIRS += examples
#if AMDISTRIBUTESOFTWARE
#SUBDIRS +=
#else
SUBDIRS += tests
#endif

if AMHAVEDOXYGEN
   SUBDIRS += doc
endif

# Do not set directories in EXTRA_DIST.
# Otherwise, the whole content of the directories will be copied
# (including .svn subfolders).
# http://www.gnu.org/software/automake/manual/html_node/Basics-of-Distribution.html
# [Last access 22.11.2013, KF]
EXTRA_DIST = LICENSE

if AMDISTRIBUTESOFTWARE
# As one cannot assume that the files are writable, the modus of the files
# will be changed, first.
dist-hook:
	chmod u+w -R $(distdir)/config/m4
	rm -rf $(distdir)/config/m4
	chmod u+w $(distdir)/configure.ac
	rm -f $(distdir)/configure.ac
	chmod u+w $(distdir)/aclocal.m4
	rm -f $(distdir)/aclocal.m4
	find $(distdir) -type f -iname \*.am -exec chmod u+w {} \;
	find $(distdir) -type f -iname \*.am -exec rm -f {} \;
endif

# Remove temporary files introduced by m4 script "ax_prefix_config_h.m4".
distclean-local:
	rm -rf $(top_builddir)/_configs.sed
	rm -rf $(top_builddir)/src/ScaFES_Config.hpp