% [masks,pind] = opSFpyrInit(sz, ht, order, twidth)
%
% Construct a steerable pyramid on matrix IM, in the Fourier domain.
% This is similar to buildSpyr, except that:
%
%    + Reconstruction is exact (within floating point errors)
%    + It can produce any number of orientation bands.
%    - Typically slower, especially for non-power-of-two sizes.
%    - Boundary-handling is circular.
%
% HEIGHT (optional) specifies the number of pyramid levels to build. Default
% is maxPyrHt(size(IM),size(FILT));
%
% The squared radial functions tile the Fourier plane, with a raised-cosine
% falloff.  Angular functions are cos(theta-k\pi/(K+1))^K, where K is
% the ORDER (one less than the number of orientation bands, default= 3).
%
% TWIDTH is the width of the transition region of the radial lowpass
% function, in octaves (default = 1, which gives a raised cosine for
% the bandpass filters).
%
% PYR is a vector containing the N pyramid subbands, ordered from fine
% to coarse.  INDICES is an Nx2 matrix containing the sizes of
% each subband.  This is compatible with the MatLab Wavelet toolbox.
% See the function STEER for a description of STEERMTX and HARMONICS.

% Eero Simoncelli, 5/97.
% See http://www.cns.nyu.edu/~eero/STEERPYR/ for more
% information about the Steerable Pyramid image decomposition.

function [masks,pind] = opSFpyrInit(dims, ht, order, twidth)

%-----------------------------------------------------------------
%% DEFAULTS:

max_ht = floor(log2(min(dims))) - 2;

if (exist('ht','var') ~= 1)
  ht = max_ht;
else
  if (ht > max_ht)
    error('Cannot build pyramid higher than %d levels.',max_ht);
  end
end

if (exist('order','var') ~= 1)
  order = 3;
elseif ((order > 15)  || (order < 0))
  fprintf(1,'Warning: ORDER must be an integer in the range [0,15]. Truncating.\n');
  order = min(max(order,0),15);
else
  order = round(order);
end
nbands = order+1;

if (exist('twidth','var') ~= 1)
  twidth = 1;
elseif (twidth <= 0)
  fprintf(1,'Warning: TWIDTH must be positive.  Setting to 1.\n');
  twidth = 1;
end

%-----------------------------------------------------------------

ctr = ceil((dims+0.5)/2);

[xramp,yramp] = meshgrid( ([1:dims(2)]-ctr(2))./(dims(2)/2), ...
    ([1:dims(1)]-ctr(1))./(dims(1)/2) );
angle = atan2(yramp,xramp);
log_rad = sqrt(xramp.^2 + yramp.^2);
log_rad(ctr(1),ctr(2)) =  log_rad(ctr(1),ctr(2)-1);
log_rad  = log2(log_rad);

%% Radial transition function (a raised cosine in log-frequency):
[Xrcos,Yrcos] = rcosFn(twidth,(-twidth/2),[0 1]);
Yrcos = sqrt(Yrcos);

YIrcos = sqrt(1.0 - Yrcos.^2);
lo0mask = pointOp(log_rad, YIrcos, Xrcos(1), Xrcos(2)-Xrcos(1), 0);

[masks,pind] = opSFpyrLevs(lo0mask, log_rad, Xrcos, Yrcos, angle, ht, nbands);

hi0mask = pointOp(log_rad, Yrcos, Xrcos(1), Xrcos(2)-Xrcos(1), 0);

pind = [size(hi0mask); pind];
masks=[{hi0mask};masks];


% [PYR, INDICES] = buildSFpyrLevs(LODFT, LOGRAD, XRCOS, YRCOS, ANGLE, HEIGHT, NBANDS)
%
% Recursive function for constructing levels of a steerable pyramid.  This
% is called by buildSFpyr, and is not usually called directly.

% Eero Simoncelli, 5/97.

function [pmasks,pind] = opSFpyrLevs(lo0mask,log_rad,Xrcos,Yrcos,angle,ht,nbands)

if (ht <= 0)
  pind = size(lo0mask);
  pmasks={lo0mask};
else

  bmasks=cell(nbands,1);
  bind = zeros(nbands,2);

%  log_rad = log_rad + 1;
  Xrcos = Xrcos - log2(2);  % shift origin of lut by 1 octave.

  lutsize = 1024;
  Xcosn = pi*[-(2*lutsize+1):(lutsize+1)]/lutsize;  % [-2*pi:pi]
  order = nbands-1;
  %% divide by sqrt(sum_(n=0)^(N-1)  cos(pi*n/N)^(2(N-1)) )
  %% Thanks to Patrick Teo for writing this out :)
  const = (2^(2*order))*(factorial(order)^2)/(nbands*factorial(2*order));
  Ycosn = sqrt(const) * (cos(Xcosn)).^order;
  himask = pointOp(log_rad, Yrcos, Xrcos(1), Xrcos(2)-Xrcos(1), 0);

  for b = 1:nbands
    anglemask = pointOp(angle, Ycosn, Xcosn(1)+pi*(b-1)/nbands, Xcosn(2)-Xcosn(1));
    bandmask = ((-sqrt(-1))^order) .* lo0mask .* anglemask .* himask;
    bmasks{b}=bandmask;
    bind(b,:)  = size(bandmask);
  end

  dims = size(lo0mask);
  ctr = ceil((dims+0.5)/2);
  lodims = ceil((dims-0.5)/2);
  loctr = ceil((lodims+0.5)/2);
  lostart = ctr-loctr+1;
  loend = lostart+lodims-1;

  log_rad = log_rad(lostart(1):loend(1),lostart(2):loend(2));
  angle = angle(lostart(1):loend(1),lostart(2):loend(2));
  YIrcos = abs(sqrt(1.0 - Yrcos.^2));
  lomask = pointOp(log_rad, YIrcos, Xrcos(1), Xrcos(2)-Xrcos(1), 0);

  [nmasks,nind] = opSFpyrLevs(lomask, log_rad, Xrcos, Yrcos, angle, ht-1, nbands);

  pmasks = [bmasks; nmasks];
  pind = [bind; nind];
   
end


