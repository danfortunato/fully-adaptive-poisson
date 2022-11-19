-include make.inc

default: mwrap mex

# Rebuild constructOperators.m and gateway.cpp via MWrap
mwrap: constructOperators.m gateway.cpp
constructOperators.m: constructOperators.mw
	$(MWRAP) -mex gateway -m constructOperators.m constructOperators.mw
gateway.cpp: constructOperators.mw
	$(MWRAP) -mex gateway -c gateway.cpp constructOperators.mw

# Rebuild gateway.mex* via MEX
mex: gateway.mex*
gateway.mex*: gateway.cpp
	$(MEX) gateway.cpp $(INCDIR) $(MEXFLAGS) $(LIBDIR) $(LIBS) -output gateway

# Remove the MEX interface
clean:
	rm -f gateway.mex*

# Remove the MEX interface, MATLAB caller, and MEX file
# Note: You will need MWrap to rebuild the deleted files!
mwrapclean: clean
	rm -f constructOperators.m gateway.cpp

.PHONY: mwrap mex clean mwrapclean
