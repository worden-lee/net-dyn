VXLDIR = ../vxl
ESTRDIR = ../libexecstream
# add -pg for profiling
CFLAGS=-g -I$(VXLDIR)/core -I$(VXLDIR)/vcl -I$(VXLDIR)/core/vnl -I$(VXLDIR)/v3p/netlib -I$(VXLDIR)/v4p/netlib -I$(ESTRDIR)
LDFLAGS=$(ESTRDIR)/exec-stream.o -L$(VXLDIR)/core/vnl/algo -L$(VXLDIR)/core/vnl/ -L$(VXLDIR)/vcl -L$(VXLDIR)/v3p/netlib -L$(VXLDIR)/lib -lvnl_algo -lvnl -lvcl -lnetlib -lv3p_netlib -lpthread

# enable to use compiler optimization
OPTIMIZE=yes
ifeq ($(OPTIMIZE),yes)
CFLAGS+=-O3
endif

#default: swarming
default: network-optimize

over: clean default

NOSRC = network-optimize.cpp \
	TimeSeriesDisplay.cpp GnuplotDisplay.cpp Parameters.cpp
NOOBJ = $(NOSRC:.cpp=.o)

SWSRC = swarming.cpp  \
	TimeSeriesDisplay.cpp GnuplotDisplay.cpp Parameters.cpp
SWOBJ = $(SWSRC:.cpp=.o)

analyse-net : analyse-net.c
	$(CC) $(CFLAGS) -ligraph $< -o $@ 

%/betweenness.csv : %/dot/1.edgelist ./analyse-net
	./analyse-net $< $@

%.edgelist : %.dot ./dot2edgelist
	./dot2edgelist $*.dot -o $@

# %.o : %.cpp
# 	$(CXX) $(CFLAGS) $(CXXFLAGS) $< -c -o $@

# fancy GNU-style line for tracking header dependencies in .P files
%.o : %.cpp
	$(CXX) $(CFLAGS) $(CXXFLAGS) -MD -c $< -o $@
	@cp $*.d $*.P; \
	    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	        -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	    rm -f $*.d

NODEP = $(NOOBJ:.o=.P)

-include $(NODEP)

SWDEP = $(SWOBJ:.o=.P)

-include $(SWDEP)

network-optimize : $(NOOBJ)
	$(CXX) $(CFLAGS) $(CXXFLAGS) -o $@ $(NOOBJ) $(LDFLAGS)

swarming : $(SWOBJ)
	$(CXX) $(CFLAGS) $(CXXFLAGS) -o $@ $(SWOBJ) $(LDFLAGS)

clean:
	$(RM) *.o *.P *~ *#

clear:
	$(RM) -r out

%.eps : %.dot
	neato -Goverlap=scale -Gsplines=true -Tps -o $@ $<

# when network.dot is changed to network~2~
%~.eps : %~
	neato -Goverlap=scale -Gsplines=true -Tps -o $@ $<
