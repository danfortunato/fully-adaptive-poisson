-include make.inc

default: mwrap mex

# Rebuild getElemCoordinates_.m and gateway.cpp via MWrap
mwrap: getElemCoordinates_.m gateway.cpp
getElemCoordinates_.m: getElemCoordinates.mw
	$(MWRAP) -mex util.gateway -m getElemCoordinates_.m getElemCoordinates.mw
gateway.cpp: getElemCoordinates.mw
	$(MWRAP) -mex util.gateway -c gateway.cpp getElemCoordinates.mw

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
	rm -f getElemCoordinates_.m gateway.cpp

.PHONY: mwrap mex clean mwrapclean
