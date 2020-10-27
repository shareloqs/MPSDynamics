try
    MATDIR=ENV["MATDIR"]
catch
    MATDIR="~/"
end
try
    DEFSAVEDIR=ENV["DEFSAVEDIR"]
catch
    DEFSAVEDIR="~/"
end
DEFCONVTHRESH = 10^-3
DEFLCTHRESH = 10^-3 # smaller threshold means faster lightcone
DEFPREC = 10^-5
DEFDLIM = 100


