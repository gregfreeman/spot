% [PYR, INDICES, STEERMTX, HARMONICS] = buildSFpyr(IM, HEIGHT, ORDER, TWIDTH)
%
% Construct a steerable pyramid on matrix IM, in the Fourier domain.
% This is similar to buildSpyr, except that:
%

function [pyr] = opSFpyrBuild(im, levels,orientations,pind,masks)


dims = size(im);
ctr = ceil((dims+0.5)/2);

m=sum(pind(:,1).*pind(:,2));

pyr=zeros(m,1);
iBand=1;
iIdx=1;
imdft = fftshift(fft2(im));
hi0mask = masks{iBand};
hi0dft =  imdft .* hi0mask;
hi0 = ifft2(ifftshift(hi0dft));
n=numel(hi0dft);
pyr(iIdx:iIdx+n-1) = real(hi0(:));
iIdx=iIdx+n;
iBand=iBand+1;

opSFpyrLevsBuild(imdft, levels, orientations,iBand,iIdx);



% [PYR, INDICES] = buildSFpyrLevs(LODFT, LOGRAD, XRCOS, YRCOS, ANGLE, HEIGHT, NBANDS)
%
% Recursive function for constructing levels of a steerable pyramid.  This
% is called by buildSFpyr, and is not usually called directly.

% Eero Simoncelli, 5/97.

    function [iBand,iIdx] = opSFpyrLevsBuild(lodft,ht,nbands,iBand,iIdx)

        n=numel(lodft);
        if (ht <= 0)
            lo0mask=masks{iBand};
            lo0 = ifft2(ifftshift(lodft.*lo0mask));      
            pyr(iIdx:iIdx+n-1) = real(lo0(:));
            iIdx=iIdx+n;
            iBand=iBand+1;

        else

          for b = 1:nbands
            banddft = lodft .* masks{iBand};
            band = ifft2(ifftshift(banddft));

            pyr(iIdx:iIdx+n-1) = real(band(:));
            iIdx=iIdx+n;
            iBand=iBand+1;
          end

          dims = size(lodft);
          ctr = ceil((dims+0.5)/2);
          lodims = ceil((dims-0.5)/2);
          loctr = ceil((lodims+0.5)/2);
          lostart = ctr-loctr+1;
          loend = lostart+lodims-1;

          lodft = lodft(lostart(1):loend(1),lostart(2):loend(2));


          [iBand,iIdx] = opSFpyrLevsBuild(lodft,ht-1,nbands,iBand,iIdx);


        end
    end
end
