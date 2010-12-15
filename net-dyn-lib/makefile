NETDYNSRC = Parameters.cpp GnuplotDisplay.cpp TimeSeriesDisplay.cpp
NETDYNOBJ = $(NETDYNSRC:.cpp=.o)
NETDYNLIB = libnet-dyn.a

default: $(NETDYNLIB)

$(NETDYNLIB) : $(NETDYNOBJ)
	$(AR) rs $(NETDYNLIB) $(NETDYNOBJ)

VXLDIR ?= ../vxl
ESTRDIR ?= ../libexecstream
# add -pg for profiling
CFLAGS=-g -I$(NETDYNDIR) -I$(VXLDIR)/core -I$(VXLDIR)/vcl -I$(VXLDIR)/core/vnl -I$(VXLDIR)/v3p/netlib -I$(VXLDIR)/v4p/netlib -I$(ESTRDIR)
LDFLAGS=-L$(NETDYNDIR) -L$(VXLDIR)/core/vnl/algo -L$(VXLDIR)/core/vnl/ -L$(VXLDIR)/vcl -L$(VXLDIR)/v3p/netlib -L$(VXLDIR)/lib $(ESTRDIR)/exec-stream.o -lnet-dyn -lvnl_algo -lvnl -lvcl -lnetlib -lv3p_netlib -lpthread

# stuff for automatic header file dependencies
%.o : %.cpp
	$(CXX) $(CFLAGS) -MD -c $< -o $@
	@cp $*.d $*.P; \
	    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	        -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	    rm -f $*.d

%.o : %.c
	$(CC) $(CFLAGS) -MD -c $< -o $@
	@cp $*.d $*.P; \
	    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	        -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	    rm -f $*.d

%.E : %.cpp
	$(CXX) $(CFLAGS) -E $< -o $@

PXXS = $(NETDYNSRC:.cpp=.P)
PFILES = $(PXXS:.c=.P)
-include $(PFILES)
