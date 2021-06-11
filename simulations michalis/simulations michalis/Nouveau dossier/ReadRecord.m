function [RecData,vel,disp]=ReadRecord(fname,recdir,rectype)
%
% Read records
% [RecData,vel,disp]=ReadRecords(fname,recdir,rectype)
% defaults: recdir = current dir, rectype = 'NGA'
% recType available = 'Silva', 'NGA', 'Somerville', 'dominic', 'Vicky','ic'
%
% RecData.rec in g
% RecData.nsteps
% RecData.dtacc
% vel (optional) in m/s
% disp (optional) in m

if isempty(recdir)
    recdir=cd;
end;
curdir=cd;
cd(recdir);

% READ THE RECORD
fid=fopen(fname);

if strcmp(rectype,'NGA')
    
    fid=fopen(fname);
    
    if fid<0
        error('##### unable to open record file ####')
    end;
    
    skipline=3;
    nsteps=0;
    dtacc=0;
    for i=1:skipline+1
        tline = fgetl(fid);
    end;
    match1 = find(tline==' ');
    match2 = find(tline=='N');
    
    nsteps=str2num(tline(1:match1(1)));
    dtacc=str2num(tline(match1(1)+1:match2-1));
    
    if nsteps<=0
        disp('ERROR: not able to find DT or NPTS');
        return;
    end;
    if dtacc<=0
        disp('ERROR: not able to find DT or NPTS');
        return;
    end;
    
    frewind(fid);
    icount=0;
    idt=0;
    acc=zeros(nsteps,1);
    while feof(fid) == 0
        tline = fgetl(fid);
        icount=icount+1;
        if icount > skipline+1
            
            % do the search in the string
            aa=isspace(tline);
            daa=zeros(length(aa)-1,1);
            for i=1:(length(aa)-1)
                daa(i)=aa(i+1)-aa(i);
            end
            daa=[daa ; 1];
            % to sting arxizei sto -1 kai teleiwnei sto 1
            ss=find(daa==-1);
            ee=find(daa==1);
            
            for i=1:length(ss)
                idt=idt+1;
                acc(idt)=str2double(tline(ss(i):ee(i)));
            end;
            
        end;
    end;
    
elseif strcmp(rectype,'NGA2')
    
    %   PEER NGA Rotated Accelerogram (v1.1)
    %   H1 for rotation: IMPERIAL VALLEY 10/15/79 2316, EC CO CENTER FF, 092 (CDMG STATION 5154)
    %   H2 for rotation: IMPERIAL VALLEY 10/15/79 2316, EC CO CENTER FF, 002 (CDMG STATION 5154)
    %   rotation angle - clockwise  51.000000
    %   SN component, azimuth =  53.000000
    %   Accelerations in g
    %   7997 0.005000 NPTS, DT
    %   0.000000 	 0.001783 	 0.001777 	 0.001760 	 0.001744
    %   0.001725 	 0.001718 	 0.001709 	 0.001707 	 0.001701
    
    fid=fopen(fname);
    
    if fid<0
        error('##### unable to open record file ####')
    end;
    
    skipline=6;
    nsteps=0;
    dtacc=0;
    for i=1:skipline+1
        tline = fgetl(fid);
    end;
    match1 = find(tline==' ');
    match2 = find(tline=='N');
    
    nsteps=str2num(tline(1:match1(1)));
    dtacc=str2num(tline(match1(1)+1:match2-1));
    
    if nsteps<=0
        disp('ERROR: not able to find DT or NPTS');
        return;
    end;
    if dtacc<=0
        disp('ERROR: not able to find DT or NPTS');
        return;
    end;
    
    frewind(fid);
    icount=0;
    idt=0;
    acc=zeros(nsteps,1);
    while feof(fid) == 0
        tline = fgetl(fid);
        icount=icount+1;
        if icount > skipline+1
            
            % do the search in the string
            aa=isspace(tline);
            daa=zeros(length(aa)-1,1);
            for i=1:(length(aa)-1)
                daa(i)=aa(i+1)-aa(i);
            end
            daa=[daa ; 1];
            % to sting arxizei sto -1 kai teleiwnei sto 1
            ss=find(daa==-1);
            ee=find(daa==1);
            
            for i=1:length(ss)
                idt=idt+1;
                acc(idt)=str2double(tline(ss(i):ee(i)));
            end;
            
        end;
    end;
    
    
elseif strcmp(rectype,'Silva')
    
    %     0.01
    %     PACIFIC ENGINEERING AND ANALYSIS STRONG-MOTION DATA
    %     IMPERIAL VALLEY 10/15/79 2316, COMPUERTAS, 285
    %     ACCELERATION TIME HISTORY IN UNITS OF G
    %     NPTS=  3600, DT= .01000 SEC
    %     -.1051707E-04   .3752166E-05   .2955270E-05  -.9803295E-04  -.2025800E-03
    
    % READ THE RECORD
    fid=fopen(fname);
    skipline=0;
    nsteps=0;
    dtacc=0;
    while feof(fid) == 0
        tline = fgetl(fid);
        match1 = findstr(tline,'NPTS=');
        num1 = length(match1);
        skipline=skipline+1;
        if num1 > 0
            nsteps=str2double(tline(6:11));
            dtacc=str2double(tline(17:24));
            break;
        end;
    end;
    
    if nsteps<=0
        disp('ERROR: not able to find DT or NPTS');
        return;
    end;
    if dtacc<=0
        disp('ERROR: not able to find DT or NPTS');
        return;
    end;
    
    frewind(fid);
    icount=0;
    idt=0;
    acc=zeros(nsteps,1);
    while feof(fid) == 0
        tline = fgetl(fid);
        icount=icount+1;
        if icount > skipline
            
            % do the search in the string
            aa=find(tline==' ');
            aa=[aa length(tline)];
            
            for i=1:length(aa)-1
                dd=tline(aa(i):aa(i+1));
                if length(dd) > 2
                    idt=idt+1;
                    acc(idt)=str2double(dd);
                end;
            end;
        end;
    end;
    
    
    
    %%%%%%%%%%%% 'European Ground Motion Database'  %%%%%%%%%%
elseif strcmp(rectype,'ic')
    
    %     sampling period:               0.010000s
    %     number of samples:             2499
    %     record length:                 24.980s
    %     units:                        m/s*s, m/s, m & s
    
    fid=fopen(fname);
    skipline=0;
    nsteps=0;
    dtacc=0;
    while feof(fid) == 0
        tline = fgetl(fid);
        match1 = findstr(tline,'sampling period:');
        num1 = length(match1);
        match2 = findstr(tline,'number of samples:');
        num2 = length(match2);
        match3 = findstr(tline,'-> corrected acceleration time histories');
        num3 = length(match3);
        
        skipline=skipline+1;
        if num1 > 0
            i1=find(tline==':')+1';
            i2=find(tline=='s')-1';
            dtacc=str2double(tline(i1:i2(2)));
        end;
        if num2 > 0
            i1=find(tline==':')+1';
            nsteps=str2double(tline(i1:end));
        end;
        if num3 > 0
            break;
        end;
        
    end;
    
    if nsteps<=0 || dtacc<=0
        disp('ERROR: not able to find DT or NPTS');
        return;
    end;
    
    frewind(fid);
    icount=0;
    idt=0;
    acc=zeros(nsteps,1);
    while feof(fid) == 0
        tline = fgetl(fid);
        icount=icount+1;
        if icount > skipline
            
            % do the search in the string
            aa=isspace(tline);
            daa=zeros(length(aa)-1,1);
            for i=1:(length(aa)-1)
                daa(i)=aa(i+1)-aa(i);
            end
            daa=[daa ; 1];
            % to sting arxizei sto -1 kai teleiwnei sto 1
            ss=find(daa==-1);
            ee=find(daa==1);
            
            for i=1:length(ss)
                idt=idt+1;
                acc(idt)=str2double(tline(ss(i):ee(i)));
            end;
            
        end;
    end;
    
    %nsteps=length(acc);
    acc=acc/9.81;
    
    
elseif strcmp(rectype,'plain1')
    
    skipline=0;
    nsteps=0;
    dtacc=0;
    while feof(fid) == 0
        tline = fgetl(fid);
        match1 = findstr(tline,'NPTS=');
        num1 = length(match1);
        skipline=skipline+1;
        if num1 > 0
            nsteps=str2double(tline(6:11));
            dtacc=str2double(tline(17:24));
            break;
        end;
    end;
    
    if nsteps<=0
        disp('ERROR: not able to find DT or NPTS');
        return;
    end;
    if dtacc<=0
        disp('ERROR: not able to find DT or NPTS');
        return;
    end;
    
    frewind(fid);
    icount=0;
    idt=0;
    while feof(fid) == 0
        tline = fgetl(fid);
        icount=icount+1;
        if icount > skipline
            
            % do the search in the string
            aa=find(tline==' ');
            aa=[aa length(tline)];
            
            for i=1:length(aa)-1
                dd=tline(aa(i):aa(i+1));
                if length(dd) > 2
                    idt=idt+1;
                    acc(idt)=str2double(dd);
                end;
            end;
            
        end;
    end;
    
    %%%%%%%%%%%% 'plain 2'  %%%%%%%%%%%%%%%%%%%
elseif strcmp(rectype,'plain2')
    
    % READ THE RECORD
    aa=load(fname);
    nsteps=length(aa)-1;
    acc=aa(2:nsteps);
    dtacc=aa(1);
    
    
    
elseif strcmp(rectype,'plain3')
    
    fid=fopen(fname);
    if fid<0
        fname
        error('##### unable to open record file ####')
    end;
    
    tline = fgetl(fid);
    tt=find(tline=='=');
    dtacc=str2num(tline(tt+1:end));
    
    icount=0;
    while feof(fid) == 0
        tline = fgetl(fid);
        icount=icount+1;
        acc(icount)=str2double(tline);
    end;
    nsteps=icount;
    acc=acc/9.81;
    
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integratre to get velocity and displacement
% (assume that acc is always in g unints).
%
vel=cumtrapz((1:nsteps)*dtacc,acc(1:nsteps)*9.81);
disp=cumtrapz((1:nsteps)*dtacc,vel(1:nsteps));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose(fid);
cd(curdir);

RecData.rec=acc;
RecData.nsteps=nsteps;
RecData.dtacc=dtacc;
