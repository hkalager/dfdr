function [FDRreal, portsize] = ComputeRealFDR(PORT, nboutperf)

portsize = sum(PORT);
if portsize > 0
    FDRreal = sum(PORT(nboutperf+1:end)) / portsize;
else
    FDRreal = 0;
end
