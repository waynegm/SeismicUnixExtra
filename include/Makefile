# Makefile for ..SeismicUnixExtra/include
# Just say "make"

include $(CWPROOT)/src/Makefile.config

INCLUDE = $I/sux.h

INSTALL	: $(INCLUDE) $I
	-rm -f INSTALL
	touch $@

$(INCLUDE):	$(ITARGET)
	@-rm -f $@
	@cp $? $@
	@$(ICHMODLINE)
	@echo $(@F) installed in $I

remake	:
	@-rm -f $(INCLUDE) INSTALL
	@$(MAKE) INSTALL
