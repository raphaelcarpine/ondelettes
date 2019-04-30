function racines = racinesPolyFrac(coefficients, exposants, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

p = inputParser ;

qmaxDef = 100;
precisDef = false;

addRequired(p,'coefficients')
addRequired(p,'exposants')
addParameter(p,'qmax', qmaxDef);
addParameter(p,'precis', precisDef);

parse(p, coefficients, exposants, varargin{:});

qmax = p.Results.qmax ;
precis = p.Results.precis;

coeffs = coefficients;
exp = exposants;
n = length(coeffs);


%% choix de q, qui minimise la somme des carrés des ecarts

p = ones(1, n);
q = 1;
for k = 1:qmax
    l = round(exp*k);
    if sum((exp-l/k).^2) < sum((exp-p/q).^2)
        p = l;
        q = k;
    end
end

%%
Polynome = zeros(1, max(p)+1);

for i = 1:length(coeffs)
    Polynome(p(i)+1) = Polynome(p(i)+1) + coeffs(i);
end


r = roots(Polynome);

r = r(angle(r) > -pi/q);
r = r(angle(r) <= pi/q);

racines = r.^q;

%%
if precis
    %TODO
end


end

