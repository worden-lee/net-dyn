# if you run make in this file it has to invoke the subdirs'
# makefiles in order to do the right thing.

egt/% : /proc/uptime
	$(MAKE) -C egt $*

net-dyn-lib/% : /proc/uptime
	$(MAKE) -C net-dyn-lib $*

landscape/% : /proc/uptime
	$(MAKE) -C landscape $*

libexecdir/exec-stream.o : libexecdir/exec-stream.cpp libexecdir/exec-stream.h
	$(CXX) $(CPPFLAGS) $(CFLAGS) $(CXXFLAGS) -c $< -o $@
