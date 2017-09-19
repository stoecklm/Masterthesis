if AMHAVEDOXYGEN

docpkg =$(PACKAGE_TARNAME)-doxy-$(PACKAGE_VERSION).zip

doc_DATA =$(docpkg)

$(docpkg): doxygen.stamp
	zip -r - html > $@

doxygen.stamp: Doxyfile \
               $(top_builddir)/src/ScaFES.hpp \
               $(top_builddir)/examples/HeatEqnFDM/HeatEqnFDM.hpp \
               $(top_builddir)/examples/WaveEqnFDM/WaveEqnFDM.hpp \
               $(top_builddir)/examples/ZhangExampleFDM/ZhangExampleFDM.hpp \
               $(top_builddir)/examples/LameNavierEqnFDM/LameNavierEqnFDM.hpp \
               $(top_builddir)/examples/LameNavierEqnFDMFirstOrder/LameNavierEqnFDMFirstOrder.hpp
	doxygen $<
	echo Timestamp > $@

install-data-hook:
	cd $(DESTDIR)$(docdir) && unzip $(docpkg)

uninstall-hook:
	cd $(DESTDIR)$(docdir) && rm -rf html

CLEANFILES = doxywarn.txt doxygen.stamp $(docpkg)

clean-local:
	rm -rf html
	rm -rf latex
endif
