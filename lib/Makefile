# Makefile for ...SeismicUnixExtra/lib

include $(CWPROOT)/src/Makefile.config

D = $I/sux.h

LIB = $L/libsux.a

ARCH =			\
	$(LIB)(sdft.o)	\
	$(LIB)(sdct.o)	\
	$(LIB)(otrcbuf.o) \
	$(LIB)(ctrcbuf.o)

INSTALL:	$(LIB) $L
	@-rm -f INSTALL
	@touch $@

$(ARCH)	: $D

$(LIB)	:	$(ARCH)
	$(RANLIB) $(LIB)

remake	:
	@-rm -f $(LIB) INSTALL
	@$(MAKE) INSTALL

list	:
	$(AR) tv $(LIB)

clean::
	rm -f junk* core a.out

.PRECIOUS:	 $(LIB)
