function poles = getPoles3ddl(f1, f2, mu1, mu2, fTMD, zetaTMD)
%GETPOLES3DDL Summary of this function goes here
%   Detailed explanation goes here

w1 = 2*pi*f1;
w2 = 2*pi*f2;
w0 = 2*pi*fTMD;
z = zetaTMD;


detM = [ w0^2*w1^2*w2^2, 2*w0*w1^2*w2^2*z, w0^2*w1^2 + w0^2*w2^2 + w1^2*w2^2 + mu1*w0^2*w2^2 + mu2*w0^2*w1^2, 2*w0*w1^2*z + 2*w0*w2^2*z + 2*mu1*w0*w2^2*z + 2*mu2*w0*w1^2*z, mu1*w0^2 + mu2*w0^2 + w0^2 + w1^2 + w2^2, 2*w0*z + 2*mu1*w0*z + 2*mu2*w0*z, 1];
detM = fliplr(detM);

poles = roots(detM);

poles = sort(poles);
poles = transpose(poles);

end



% syms p w1 w2 mu1 mu2 w0 z
% 
% M = [p^2 + w1^2, 0, mu1*p^2;
%     0, p^2 + w2^2, mu2*p^2;
%     -w0^2 - 2*z*w0*p, -w0^2 - 2*z*w0*p, p^2 + w0^2 + 2*z*w0*p];
% 
% d = det(M)
% 
% coeffs(d, p)
% 
% [ w0^2*w1^2*w2^2, 2*w0*w1^2*w2^2*z, w0^2*w1^2 + w0^2*w2^2 + w1^2*w2^2 + mu1*w0^2*w2^2 + mu2*w0^2*w1^2, 2*w0*w1^2*z + 2*w0*w2^2*z + 2*mu1*w0*w2^2*z + 2*mu2*w0*w1^2*z, mu1*w0^2 + mu2*w0^2 + w0^2 + w1^2 + w2^2, 2*w0*z + 2*mu1*w0*z + 2*mu2*w0*z, 1]
%  