% RES = reconSFpyr(PYR, INDICES, LEVS, BANDS, TWIDTH)
%
% Reconstruct image from its steerable pyramid representation, in the Fourier
% domain, as created by buildSFpyr.
%
% PYR is a vector containing the N pyramid subbands, ordered from fine
% to coarse.  INDICES is an Nx2 matrix containing the sizes of
% each subband.  This is compatible with the MatLab Wavelet toolbox.
%
% LEVS (optional) should be a list of levels to include, or the string
% 'all' (default).  0 corresonds to the residual highpass subband.  
% 1 corresponds to the finest oriented scale.  The lowpass band
% corresponds to number spyrHt(INDICES)+1.
%
% BANDS (optional) should be a list of bands to include, or the string
% 'all' (default).  1 = vertical, rest proceeding anti-clockwise.
%
% TWIDTH is the width of the transition region of the radial lowpass
% function, in octaves (default = 1, which gives a raised cosine for
% the bandpass filters).

%%% MODIFIED VERSION, 7/04, uses different lookup table for radial frequency!

% Eero Simoncelli, 5/97.

function res = opSFpyrRecon(pyr, levels,orientations,pind,masks)

dims = pind(1,:);
ctr = ceil((dims+0.5)/2);
iBand=1;
iIdx=1;


% residual highpass subband
sz=pind(iBand,:);
n=prod(sz);
hidft = fftshift(fft2(reshape(pyr(iIdx:iIdx+n-1),sz(1),sz(2))));
resdft = hidft .*  masks{iBand};
iIdx=iIdx+n;
iBand=iBand+1;
% recurse bands

opSFpyrLevsRecon(levels);
% ifft

res = real(ifft2(ifftshift(resdft)));

    % RESDFT = reconSFpyrLevs(PYR,INDICES,LOGRAD,XRCOS,YRCOS,ANGLE,NBANDS,LEVS,BANDS)
    %
    % Recursive function for reconstructing levels of a steerable pyramid
    % representation.  This is called by reconSFpyr, and is not usually
    % called directly.

    % Eero Simoncelli, 5/97.

    function opSFpyrLevsRecon(ht)

        sz = pind(iBand,:);
        n=prod(sz);
        loctr = ceil((sz+0.5)/2);
        lostart = ctr-loctr+1;
        loend = lostart+sz-1;

        if (ht <= 0)
          nresdft = fftshift(fft2(reshape(pyr(iIdx:iIdx+n-1),sz(1),sz(2))));
          resdft(lostart(1):loend(1),lostart(2):loend(2)) = ...
                resdft(lostart(1):loend(1),lostart(2):loend(2)) + nresdft .* masks{iBand};
%           pdft = nresdft .* masks{iBand};
          iIdx=iIdx+n;
          iBand=iBand+1;
        else
          pdft=zeros(sz);
          for b = 1:orientations
              nresdft = fftshift(fft2(reshape(pyr(iIdx:iIdx+n-1),sz(1),sz(2))));

%              resdft(lostart(1):loend(1),lostart(2):loend(2)) = ...
%                 resdft(lostart(1):loend(1),lostart(2):loend(2)) + nresdft .* masks{iBand};
              pdft = pdft-nresdft .* masks{iBand};
              iIdx=iIdx+n;
              iBand=iBand+1;
          end

          opSFpyrLevsRecon(ht-1);
            resdft(lostart(1):loend(1),lostart(2):loend(2)) = ...
                resdft(lostart(1):loend(1),lostart(2):loend(2)) + pdft;
        end

%         imagesc(real(ifft2(ifftshift(resdft))))
    end
end
