function output=bound(input,mini,maxi)
output = min(max(input,mini),maxi);
% On limite input a rester entre mini et maxi
% input = mini si input<mini
% input = maxi si input>maxi