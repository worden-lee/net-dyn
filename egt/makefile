VXLDIR = ../vxl
CFLAGS=-O3 -g -I$(VXLDIR)/core -I$(VXLDIR)/vcl -I$(VXLDIR)/core/vnl
LDFLAGS=-L$(VXLDIR)/core/vnl/algo -L$(VXLDIR)/core/vnl/ -L$(VXLDIR)/vcl -L$(VXLDIR)/v3p/netlib -lvnl_algo -lvnl -lvcl -lnetlib -lv3p_netlib

default: network-optimize

NOSRC = network-optimize.cpp DotDisplay.cpp
NOOBJ = $(NOSRC:.cpp=.o)

# %.o : %.cpp
# 	$(CXX) $(CFLAGS) $(CXXFLAGS) $< -c -o $@

# fancy GNU-style line for tracking header dependencies
%.o : %.cpp
	$(CXX) $(CFLAGS) $(CXXFLAGS) -MD -c $< -o $@
	@cp $*.d $*.P; \
	    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	        -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	    rm -f $*.d

NODEP = $(NOSRC:.cpp=.P)

-include $(NODEP)

network-optimize : $(NOOBJ)
	$(CXX) $(CFLAGS) $(CXXFLAGS) -o $@ $(NOOBJ) $(LDFLAGS)

clean:
	$(RM) *.o *~ *#

%.eps : %.dot
	neato -Goverlap=scale -Gsplines=true -Tps -o $@ $<
