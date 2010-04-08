NETDYNSRC = Parameters.cpp GnuplotDisplay.cpp TimeSeriesDisplay.cpp
NETDYNOBJ = $(NETDYNSRC:.cpp=.o)
NETDYNLIB = libnet-dyn.a

default: $(NETDYNLIB)

$(NETDYNLIB) : $(NETDYNOBJ)
	$(AR) rs $(NETDYNLIB) $(NETDYNOBJ)

CFLAGS=-g -I../libexecstream -I../vxl/core -I../vxl/vcl -I../vxl/core/vnl -I../vxl/v3p/netlib -I../vxl/v4p/netlib
LDFLAGS=../libexecstream/exec-stream.o -L../vxl/core/vnl/algo -L../vxl/core/vnl -L../vxl/vcl -L../vxl/v3p/netlib -L../vxl/lib -lnet-dyn -lvnl_algo -lvnl -lvcl -lnetlib -lv3p_netlib -lpthread

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
