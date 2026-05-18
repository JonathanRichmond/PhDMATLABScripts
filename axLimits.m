%%% axLimits
%%% Jonathan LeFevre Richmond
%%% C: 20 April 2026

function lim = axLimits(varargin)

pad = varargin{end};
allData = vertcat(varargin{1:end-1});
allData = allData(isfinite(allData));
mn = min(allData);
mx = max(allData);
r = mx-mn;
if r == 0
    r = 1;
end

lim = [mn-pad*r, mx+pad*r];