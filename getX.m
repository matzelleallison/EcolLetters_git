% get food value for desired f value
% f = desired functional response
% K = half sat coeff
function X = getX(f,K)
X = f*K/(1-f);
