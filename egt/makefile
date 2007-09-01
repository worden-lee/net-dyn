CFLAGS=-g

default: network-optimize

NOOBJ = network-optimize.o DotDisplay.o

%.o : %.cpp
	$(CXX) $(CFLAGS) $(CXXFLAGS) $< -c -o $@

network-optimize : $(NOOBJ)
	$(CXX) $(CFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $(NOOBJ)

clean:
	$(RM) *.o *~ *#
