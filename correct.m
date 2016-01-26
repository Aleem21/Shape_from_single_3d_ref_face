function corrected = correct(I,correction)
sz = size(I);
corrected = reshape((correction*(reshape(I,[],3)'))',sz);
corrected = max(0,min(corrected,1));
corrected(isnan(I)) = nan;
end