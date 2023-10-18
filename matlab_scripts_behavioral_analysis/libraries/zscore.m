function [out]= zscore_dF(in)
    mu=mean(in(1:10));
    sigma=std(in(1:10));
    out=(in-mu)/sigma;
end

