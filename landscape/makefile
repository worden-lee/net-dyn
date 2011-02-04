VXLDIR = ../vxl
ESTRDIR = ../libexecstream
# add -pg for profiling
#CFLAGS=-g -I$(VXLDIR)/core -I$(VXLDIR)/vcl -I$(VXLDIR)/core/vnl -I$(VXLDIR)/v3p/netlib -I$(VXLDIR)/v4p/netlib -I$(ESTRDIR)
#LDFLAGS=$(ESTRDIR)/exec-stream.o -L$(VXLDIR)/core/vnl/algo -L$(VXLDIR)/core/vnl/ -L$(VXLDIR)/vcl -L$(VXLDIR)/v3p/netlib -L$(VXLDIR)/lib -lvnl_algo -lvnl -lvcl -lnetlib -lv3p_netlib -lpthread
CFLAGS=-g
LDFLAGS=-lssl

# enable to use compiler optimization
OPTIMIZE=yes
ifeq ($(OPTIMIZE),yes)
CFLAGS+=-O3
endif

default: search-landscape

over: clean default

SLSRC = search-landscape.cpp  \
	BlockFitnessLandscape.cpp Point.cpp ../network/Parameters.cpp
SLOBJ = $(SLSRC:.cpp=.o)

# %.o : %.cpp
# 	$(CXX) $(CFLAGS) $(CXXFLAGS) $< -c -o $@

# fancy GNU-style line for tracking header dependencies in .P files
%.o : %.cpp
	$(CXX) $(CFLAGS) $(CXXFLAGS) -MD -c $< -o $@
	@cp $*.d $*.P; \
	    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	        -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	    rm -f $*.d

../network/% : /proc/uptime
	$(MAKE) -C ../network $*

SLDEP = $(SLOBJ:.o=.P)

-include $(SLDEP)

search-landscape : $(SLOBJ)
	$(CXX) $(CFLAGS) $(CXXFLAGS) -o $@ $(SLOBJ) $(LDFLAGS)

clean:
	$(RM) *.o *.P *~ *#

clear:
	$(RM) -r out

%.eps : %.dot
	neato -Goverlap=scale -Gsplines=true -Tps -o $@ $<
