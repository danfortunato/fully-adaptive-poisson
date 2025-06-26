-include make.inc

default: mwrap mex

# Rebuild getElemCoordinates_.m and getElemCoordinates_.cpp via MWrap
mwrap: getElemCoordinates_.m getElemCoordinates_.cpp
getElemCoordinates_.m: getElemCoordinates_.mw
	$(MWRAP) -mex util.getElemCoordinates_mex -m getElemCoordinates_.m getElemCoordinates_.mw
getElemCoordinates_.cpp: getElemCoordinates_.mw
	$(MWRAP) -mex util.getElemCoordinates_mex -c getElemCoordinates_.cpp getElemCoordinates_.mw

# Rebuild getElemCoordinates_mex.mex* via MEX
mex: getElemCoordinates_mex.mex*
getElemCoordinates_mex.mex*: getElemCoordinates_.cpp
	$(MEX) getElemCoordinates_.cpp $(INCDIR) $(MEXFLAGS) $(LIBDIR) $(LIBS) -output getElemCoordinates_mex

# Remove the MEX interface
clean:
	rm -f getElemCoordinates_mex.mex*

# Remove the MEX interface, MATLAB caller, and MEX file
# Note: You will need MWrap to rebuild the deleted files!
mwrapclean: clean
	rm -f getElemCoordinates_.m getElemCoordinates_.cpp

.PHONY: mwrap mex clean mwrapclean
