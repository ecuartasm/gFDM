if ispc
    mex BioCalc.cpp
elseif ismac    
    mex BioCalcUnix.cpp
elseif isunix
    mex BioCalcUnix.cpp
end