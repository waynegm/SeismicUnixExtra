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
	$B/suctrcmedian \
	$B/susdft_denoise \
	$B/sulpasmooth \
	$B/suvpef

INSTALL	:	$(PROGS)
	@-rm -f INSTALL
	@cat ../Stub.md ../programs.md >../README.md 
	@cat ../docs/Stub.md ../docs/programs.md > ../docs/start.md
	@-rm ../programs.md
	@-rm ../docs/programs.md
	@touch $@


$(PROGS):	$(CTARGET) $D 
	-$(CC) $(CFLAGS) $(@F).c $(LFLAGS) -o $@
	@$(MCHMODLINE)
	@-$(@F) 2> ../docs/$(@F).md
	@echo \| [$(@F)]"(docs/"$(@F)".md)" \| $(shell $(@F) 2>&1|head -n 2|tail -n 1) \| >> ../programs.md
	@echo \| [$(@F)]"("$(@F)".md)" \| $(shell $(@F) 2>&1|head -n 2|tail -n 1) \| >> ../docs/programs.md
	@echo $(@F) installed in $B

remake	:
	-rm -f $(PROGS) INSTALL
	$(MAKE) 
	
clean::
	rm -f a.out junk* JUNK* core
