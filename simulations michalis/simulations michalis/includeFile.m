%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INITIALIZE AND READ THE GROUND MOTION RECORD
%
% created by Michalis Fragiadakis, Dec 2013
% mfrag@mail.ntua.gr
%
% please report any bugs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;


% READ GROUND MOTION OF DIFFERENT FORMATS
icase=5;

recsPath = fullfile(fileparts(matlab.desktop.editor.getActiveFilename), 'recs');

if icase ==1
    [acc,vel]=ReadRecord('184057.AT2',recsPath,'NGA');
    
elseif icase ==2
    [acc,vel]=ReadRecord('PL_1_SN.acc',recsPath,'NGA2');
    
elseif icase ==3
    [acc,vel]=ReadRecord('NGA_no_28_C12050.AT2',recsPath,'NGA');
    
elseif icase ==4
    [acc,vel]=ReadRecord('AssisiStallone_X.txt',recsPath,'plain3');
    
elseif icase ==5
    [acc,vel]=ReadRecord('sv23.th',recsPath,'Silva');
    
end;


% UNCOMMENT TO PLOT THE RECORD
% figure()
% plot((1:acc.nsteps)*acc.dtacc,acc.rec)