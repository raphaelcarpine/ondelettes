function Pn = numberOfPartitions(N)
% numberOfPartitions: compute the number of partitions of the positive integer n
% usage: Pn = numberOfPartitions(N)
%
% Uses Macmahon's identity for the partitions function
% in terms of a recurrence relation based on the generalized
% pentagonal numbers.
%
%  http://mathworld.wolfram.com/PartitionFunctionP.html
%  
%  http://mathworld.wolfram.com/PentagonalNumber.html
% 
% Arguments: (input)
%  N  - a scalar integer or list of integers.
%       (Large values of N will result in HUGE
%        numbers, so beware. A practical limit
%        on the size of N is probably 
%
%  Pn - The number of partitions, P(N)
%       Pn will be a vpi number or numbers
%       if N was an array.
%
% Example:
%  Pn = numberOfPartitions(500)
% Pn =
%   2300165032574323995027
%
%  See also: 
%
%  Author: John D'Errico
%  e-mail: woodchips@rochester.rr.com
%  Release: 1.0
%  Release date: 5/10/09

if nargin ~= 1
  error('numberOfPartitions requires exactly one argument')
end

N = double(N);
Pn = vpi(N);
[N,Ntags] = sort(N(:));

if any(N<=0) || any(N~=round(N))
  error('N must be a positive integer')
end

Nmax = N(end);
% we need to generate ALL of the results
% up to Nmax. Do this sequentially using a
% recurrence relation. Rather than build them
% all up, the first 240 numbers in this sequence
% are:

Plist = [1 1 2 3 5 7 11 15 22 30 42 56 77 101 135 176 231 ...
 297 385 490 627 792 1002 1255 1575 1958 2436 3010 3718 ...
 4565 5604 6842 8349 10143 12310 14883 17977 21637 26015 ...
 31185 37338 44583 53174 63261 75175 89134 105558 124754 ...
 147273 173525 204226 239943 281589 329931 386155 451276 ...
 526823 614154 715220 831820 966467 1121505 1300156 1505499 ...
 1741630 2012558 2323520 2679689 3087735 3554345 4087968 ...
 4697205 5392783 6185689 7089500 8118264 9289091 10619863 ...
 12132164 13848650 15796476 18004327 20506255 23338469 ...
 26543660 30167357 34262962 38887673 44108109 49995925 ...
 56634173 64112359 72533807 82010177 92669720 104651419 ...
 118114304 133230930 150198136 169229875 190569292 214481126 ...
 241265379 271248950 304801365 342325709 384276336 431149389 ...
 483502844 541946240 607163746 679903203 761002156 851376628 ...
 952050665 1064144451 1188908248 1327710076 1482074143 ...
 1653668665 1844349560 2056148051 2291320912 2552338241 ...
 2841940500 3163127352 3519222692 3913864295 4351078600 ...
 4835271870 5371315400 5964539504 6620830889 7346629512 ...
 8149040695 9035836076 10015581680 11097645016 12292341831 ...
 13610949895 15065878135 16670689208 18440293320 20390982757 ...
 22540654445 24908858009 27517052599 30388671978 33549419497 ...
 37027355200 40853235313 45060624582 49686288421 54770336324 ...
 60356673280 66493182097 73232243759 80630964769 88751778802 ...
 97662728555 107438159466 118159068427 129913904637 142798995930 ...
 156919475295 172389800255 189334822579 207890420102 228204732751 ...
 250438925115 274768617130 301384802048 330495499613 362326859895 ...
 397125074750 435157697830 476715857290 522115831195 571701605655 ...
 625846753120 684957390936 749474411781 819876908323 896684817527 ...
 980462880430 1071823774337 1171432692373 1280011042268 1398341745571 ...
 1527273599625 1667727404093 1820701100652 1987276856363 2168627105469 ...
 2366022741845 2580840212973 2814570987591 3068829878530 3345365983698 ...
 3646072432125 3972999029388 4328363658647 4714566886083 5134205287973 ...
 5590088317495 6085253859260 6622987708040 7206841706490 7840656226137 ...
 8528581302375 9275102575355 10085065885767 10963707205259 ...
 11916681236278 12950095925895 14070545699287 15285151248481 ...
 16601598107914 18028182516671 19573856161145 21248279009367 ...
 23061871173849 25025873760111 27152408925615 29454549941750 ...
 31946390696157 34643126322519 37561133582570 40718063627362 ...
 44132934884255 47826239745920 51820051838712 56138148670947 ...
 60806135438329 65851585970275 71304185514919 77195892663512 ...
 83561103925871 90436839668817 97862933703585 105882246722733];

np0 = length(Plist); 
if Nmax > np0
  Plist = [Plist,repmat(vpi(0),1,Nmax - np0 + 1)];
end
for n = np0:Nmax
  % use Euler's recurrence
  p = 0;
  k = 1;
  flag = true;
  while (k <= n) && flag
    s = (-1)^(k+1);
    k1 = k*(3*k-1)/2;
    k2 = k*(3*k+1)/2;
    if k1 <= n
      p = p + s*Plist(n - k1 + 1);
    else
      flag = false;
    end
    if k2 <= n
      p = p + s*Plist(n - k2 + 1);
    end
    
    k = k + 1;
  end
  Plist(n+1) = p;
end

% extract the requested partition numbers
% from the list.
Pn(Ntags) = Plist(1+N);

