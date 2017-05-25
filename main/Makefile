# Makefile for ...su/main

include $(CWPROOT)/src/Makefile.config

D = $L/libcwp.a $L/libpar.a $L/libsu.a $L/libsux.a

LFLAGS= $(PRELFLAGS) -L$L -lsux -lsu -lpar -lcwp -lm $(POSTLFLAGS)


PROGS =			\
	$B/susdft	\
	$B/suisdft	\
	$B/susdct	\
	$B/suisdct	\
	$B/sutrcmedian	\
	$B/suctrcmedian


INSTALL	:	$(PROGS)
	@-rm -f INSTALL
	@touch $@


$(PROGS):	$(CTARGET) $D 
	-$(CC) $(CFLAGS) $(@F).c $(LFLAGS) -o $@
	@$(MCHMODLINE)
	@echo $(@F) installed in $B

remake	:
	-rm -f $(PROGS) INSTALL
	$(MAKE) 
	
clean::
	rm -f a.out junk* JUNK* core
