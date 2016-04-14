%%% BLIND INVERSE GAMMA CORRECTION
%%%    1. generate random signal from N(0,1) of length 'dim'
%%%    2. rescale signal into [0,1]
%%%    3. apply gamma correction (g(u) = u^gamma)
%%%    4. blindly estimate gamma by minimizing bicoherence
%%%       for each (possibly overlapping (wininc)) windowed signal (windim)
%%%          for each inverse gamma value in 'range'
%%%             apply inverse gamma 
%%%             compute bicoherence of inverse gamma signal
%%%          end
%%%          estimated gamma = inverse gamma with minimum bicoherence
%%%       end
%%%    5. final estimate = average of individual estimates
%%%
%%% RUN:
%%%     gammaest = blindgamma( 1.2, [0.5:0.1:2.0] );
%%%
%%% REFERENCES:
%%%    www.cs.dartmouth.edu/farid/publications/ip01.html
%%%    www.cs.dartmouth.edu/farid/publications/josa01.html
%%%

function[gammaest] = blindgamma( gamma, range,im )

dim    = 2^13;               % signal length
windim = 2^10;               % window length
wininc = 2^9;                % window increment
% f      = randn( 1, dim );    % original signal
im(isnan(im)) = rand(1,sum(isnan(im(:))));
% im(im==0) = rand(1,sum(im(:)==0));
% im = permute(im,[2 1 3]);
% im = max(im,0);
f = im(randperm(numel(im)));
% f = im(:);
dim = numel(f);
windim = size(im,1);
wininc = windim*10;
f      = f - min(f);         % normalize into [0,1]
f      = f / max(f);         % 
fg     = f.^gamma;           % apply gamma 

ind = [1:windim];
k = 1;
while( max(ind)<dim )        % for each window
   fw = fg(ind);              
   for g = 1 : length(range) % for each inverse gamma value
      bic = bispec( fw.^(1/range(g)) );
      B(k,g) = nanmean( abs(bic(:)) );
   end
   k = k + 1;
   ind = ind + wininc;
end

for k = 1 : size(B,1)
   [val,ind] = min(B(k,:));
   gammaest(k)  = range(ind);
end
gammaest = mean(gammaest);

return;


%%%
%%% COMPUTE BICOHERENCE
%%%
function [bic] = bispec( y )

[ly, nrecs] = size(y);
if( ly == 1 )
   y = y(:);
end

% ---------------- bispectrum parameters ----------------------

nsamp   = 64;
overlap = 32;
nadvance = nsamp - overlap;
nrecs   = fix ( (ly*nrecs - overlap) / nadvance);
nfft    = 128;
wind    = hanning(nsamp); 
wind    = wind(:);

% ---------------- bispectrum calculation ----------------------
bic	= zeros( nfft, nfft );
YY 	= zeros( nfft, 1 );
Yc12	= zeros( nfft, nfft );
YY12	= zeros( nfft, nfft );
mask	= hankel([1:nfft],[nfft,1:nfft-1] );
ind	= [1:nsamp];

for k = 1 : nrecs
        ys	= y(ind);
	ys	= (ys(:) - mean(ys)) .* wind;
	ys 	= ys(:) - mean(ys);
	Y	= fft(ys,nfft) / nsamp;
        Yc	= conj( Y );
	YY	= YY + ( Y .* Yc );
	Yc12(:) = Yc( mask );
	bic	= bic + (Y * Y.') .* Yc12;
        ind	= ind + nadvance;
end

YY      = sqrt( YY / nrecs );
YY12(:) = YY(mask);
bic     = bic / nrecs;
bic     = bic ./ sqrt( (abs(YY * YY.').^2) .* (abs(YY12).^2) );
bic     = fftshift( bic );

return;
