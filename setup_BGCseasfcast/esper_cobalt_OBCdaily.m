% Create monthly OBC fields of total alkalinity (TA) and
% Dissolved Inorganic Carbon DIC  for COBALT
% From daily OBC fields
%  
% Save 2 years to allow 1-yr seas f/casts with init months 1- 10
%
% Using ESPER codes
%
% In a bash script:
%#!/bin/bash
%
%# Define input arguments
%YS=1993
%YE=1995
%
%# Run MATLAB script and pass arguments as variable assignments
%matlab -nodisplay -nosplash -r "YR1=${YS}; YR2=${YE}; run('esper_cobalt_OBCdaily.m'); exit;"

addpath /home/Dmitry.Dukhovskoy/matlab/ESPER-main;
addpath /home/Dmitry.Dukhovskoy/matlab/ESPER-main/ESPER_LIR_Files;
addpath /home/Dmitry.Dukhovskoy/matlab/ESPER-main/ESPER_NN_Files;
addpath /home/Dmitry.Dukhovskoy/matlab/ESPER-main/SimpleCantEstimateFiles;
%addpath /home/Dmitry.Dukhovskoy/matlab/ESPER-main/private;
addpath /home/Dmitry.Dukhovskoy/matlab/MyMatlab;

% Read depths:
pthtopo = '/work/Dmitry.Dukhovskoy/NEP_input/topo_grid/';
pthobc = '/work/Dmitry.Dukhovskoy/NEP_input/spear_obc_daily/';
pthesper = '/work/Dmitry.Dukhovskoy/NEP_input/BGC_esper_seasfcast/';

dftopo = sprintf('%stopog.nc',pthtopo);
% To allow for the YR arguments passed from the command line:
if ~exist('YR1', 'var')
  YR1 = 2021;
end
if ~exist('YR2', 'var')
  YR2 = 2021;
end
Nyrs_save = 2;    % how many years save in 1 file
YR2    = YR2+(Nyrs_save-1);
MMI    = 1;    % SPEAR init month: 1, 4, 7, 10
enmb   = 1;    % SPEAR ens numb: 1,..., 10
requested_vars  = [1,2];  % Request total alk, DIC
predictor_types = [1,2];
equations = [8];
dnmb_ref = datenum(1993,1,1);  % reference date for NetCDF time var


DATA = struct;
for isegm = 1:4
  DATA(1).segm(isegm) = struct;
end

kYR=0;
for YR=YR1:YR2;
  kYR = kYR+1;
  ZMsegm = [];
  ZZsegm = [];

  % Create monthly means, however
  % Since the daily OBC fields are interpolated in time monthly
  % mean fields from SPEAR - not extra averaging is needed, just 
  % use the actual dates
  % include Jan 1 and Dec 31+1 to start / finish the 1st / last years 
  TM = [ ...
  datenum(YR, 1, 1), ...           
  datenum(YR, 1:12, 15), ...    
  datenum(YR, 12, 31)          
  ];

  JDAY  = TM-TM(1);
  nrecs = length(JDAY);
  TMM = TM;
  TMM(end) = TMM(end)+1; 

  for isegm=1:4;
    fprintf(' ====    %i  MMI=%2.2i ensmb=e%2.2i segment %i =====\n',YR, MMI, enmb, isegm);

    segstr = sprintf('_segment_%3.3i',isegm);
    dflnc = sprintf('%sOBCs_spear_daily_init%i%2.2i01_e%2.2i.nc',pthobc,YR,MMI,enmb);
    fprintf('Reading %s\n',dflnc);
    %ncid = netcdf.open(dflnc, 'NOWRITE');
 
    salname  = sprintf('so%s',segstr);
    tempname = sprintf('thetao%s',segstr);
    lonname  = sprintf('lon%s',segstr);
    latname  = sprintf('lat%s',segstr);
    dzname   = sprintf('dz_thetao%s',segstr);
    S3d = squeeze(ncread(dflnc, salname ));  % x,y,days since Jan 1
    T3d = squeeze(ncread(dflnc, tempname));
    lat = squeeze(ncread(dflnc, latname ));
    lon = squeeze(ncread(dflnc, lonname ));
    DZ  = squeeze(ncread(dflnc, dzname  )); % layers dz with time, write it to the output as is
    [dim1, dim2, dim3] = size(T3d);          % horiz axis, depth, days

    DATA(kYR).segm(isegm).lon = lon;
    DATA(kYR).segm(isegm).lat = lat;

    % Derive nx, ny dimensions
    if ~isfield(DATA(kYR).segm(isegm), 'nx_segm') || ...
       length(DATA(kYR).segm(isegm).nx_segm) == 0
      nxsegm_name = sprintf('nx%s',segstr);
      nysegm_name = sprintf('ny%s',segstr);
      nx_segm = ncread(dflnc, nxsegm_name);
      ny_segm = ncread(dflnc, nysegm_name);
      nx_len  = length(nx_segm);
      ny_len  = length(ny_segm);
      DATA(kYR).segm(isegm).nx_segm = nx_len;
      DATA(kYR).segm(isegm).ny_segm = ny_len;
      DATA(kYR).segm(isegm).nz_segm = dim2; 
    end
    nx = DATA(kYR).segm(isegm).nx_segm;
    ny = DATA(kYR).segm(isegm).ny_segm;
    nz = DATA(kYR).segm(isegm).nz_segm;

    % DZ is 3D array with time dim. 
    [d1,nzlev,ntime] = size(DZ);
    ZZsegm = zeros(nzlev+1,d1);
    ZMsegm = zeros(nzlev,d1);
    for kk=1:nzlev
      dh = squeeze(DZ(:,kk,1)).';
      ZZsegm(kk+1,:) = squeeze(ZZsegm(kk,:))-abs(dh);
      ZMsegm(kk,:) = 0.5*(ZZsegm(kk,:)+ZZsegm(kk+1,:));
    end

    nlong = max([nx,ny]);
    outp_coords = zeros(nlong * nzlev, 3);  
    irow = 1;
    for kk = 1:nzlev
      zmk = abs(ZMsegm(kk,:));
      outp_coords(irow:irow+nlong-1, :) = [lon(:), lat(:), zmk(:)];
      irow = irow + nlong;
    end

    % Need Jan 1 for the 1st year but not for the following years
    if kYR == 1
      kday1 = 1; 
    else
      kday1 = 2;
    end
    % Need Dec 31. only for the last record to close off the year
    if kYR < Nyrs_save
      kday2 = nrecs-1;
    else
      kday2 = nrecs;
    end
  
    for kday=kday1:kday2
      fprintf('  ---  segm %i jday=%i \n',isegm,JDAY(kday))
      est_dates = [YR + JDAY(kday)/JDAY(end)];
      % Extend into the next year for the last record
      if kYR == Nyrs_save
        est_dates = max([est_dates, YR+1+1/365]);
      end
        
      itime = JDAY(kday)+1; 
      S2d = squeeze(S3d(:,:,itime));
      T2d = squeeze(T3d(:,:,itime));

      pred_vars = zeros(nlong*nzlev, 2); 
      irow = 1;
      for k = 1:nzlev
        s = S2d(:, k);  % S at depth k
        t = T2d(:, k);  % T at depth k
        % Get rid off low S to avoid warnings:
        s(s <= 5.) = 5;
        pred_vars(irow:irow+nlong-1, :) = [s, t];  % Stack as [salt, temp]
        irow = irow + nlong;
      end

      % ESPER_NN looks similar to ESPER_LIR except for the near-coast regions
      % with low S, EPSR_NN does not have nans in bottom grid cells
      %ESPER = ESPER_LIR(requested_vars, outp_coords, pred_vars, ...
      %                  predictor_types, 'Equations', equations, ...
      %                  'EstDates', est_dates);

      ESPER = ESPER_NN(requested_vars, outp_coords, pred_vars, ...
                        predictor_types, 'Equations', equations, ...
                        'EstDates', est_dates);

      % Reshape and orient depth in rows:
      alk = (reshape(ESPER.TA, nlong, nzlev) * 1e-6)';
      dic = (reshape(ESPER.DIC, nlong, nzlev) * 1e-6)';

      % Forward fill missing values down the rows (fill bottom if NaNs):
      alk = fillmissing(alk, 'previous', 1);
      dic = fillmissing(dic, 'previous', 1);
      
      DATA(kYR).segm(isegm).alk(kday,:,:) = alk;
      DATA(kYR).segm(isegm).dic(kday,:,:) = dic;
      DATA(kYR).segm(isegm).jdays(kday)   = JDAY(kday);
      DATA(kYR).segm(isegm).dnmb(kday)    = TMM(kday);

    end  % for kday
  end    % for isegm
    
  if kYR == Nyrs_save
    % Save data
    yr_start=YR-Nyrs_save+1;
    yr_end=YR;
    flnm_esper = sprintf('bgc_esper_SPEARmnth_%i-%i.nc',yr_start,yr_end);
    flesper_out = fullfile(pthesper, flnm_esper);
    if ~exist(pthesper, 'dir')
      mkdir(pthesper);
    end

    if isfile(flesper_out)
      delete(flesper_out);
    end
    fprintf('Saving fields --> %s\n',flesper_out);

    % create netcdf, define dims:
    for isegm=1:4
      [ndays, zdim, hdim] = size(DATA(1).segm(isegm).alk);

      if isegm == 1 || isegm == 3
        xdim = hdim;
        ydim = 1;
      else
        xdim = 1;
        ydim = hdim;
      end

      segstr = sprintf('_segment_%3.3i',isegm);
      alk_var = sprintf('alk%s',segstr);
      dic_var = sprintf('dic%s',segstr);
      dzalk_var = sprintf('dz_alk%s',segstr);
      dzdic_var = sprintf('dz_dic%s',segstr);
      lat_var   = sprintf('lat%s',segstr);
      lon_var   = sprintf('lon%s',segstr);
      zdim_name = sprintf('nz%s',segstr);
      ydim_name = sprintf('ny%s',segstr);
      xdim_name = sprintf('nx%s',segstr);
     

      nccreate(flesper_out, alk_var, ...
       'Dimensions', {'time',Inf, zdim_name,zdim, ydim_name,ny, xdim_name,nx},...
       'Datatype', 'double',...
       'FillValue', 1e20);

      nccreate(flesper_out, dic_var, ...
       'Dimensions', {'time',Inf, zdim_name,zdim, ydim_name,ny, xdim_name,nx},...
       'Datatype', 'double',...
       'FillValue', 1e20);

      nccreate(flesper_out, dzalk_var, ...
       'Dimensions', {'time',Inf, zdim_name,zdim, ydim_name,ny, xdim_name,nx},...
       'Datatype', 'double',...
          'FillValue', 1e20);

      nccreate(flesper_out, dzdic_var, ...
       'Dimensions', {'time',Inf, zdim_name,zdim, ydim_name,ny, xdim_name,nx},...
       'Datatype', 'double',...
          'FillValue', 1e20);

      if isegm == 1 || isegm == 3
        nccreate(flesper_out, lat_var, ...
         'Dimensions', {xdim_name,nx},...
         'Datatype', 'double',...
          'FillValue', 1e20);

        nccreate(flesper_out, xdim_name,...
         'Dimensions', {xdim_name,nx},...
         'Datatype', 'int32',...
          'FillValue', int32(1e20));
      else
        nccreate(flesper_out, lon_var, ...
         'Dimensions', {ydim_name,ny},...
         'Datatype', 'double',...
          'FillValue', 1e20);

        nccreate(flesper_out, ydim_name,...
         'Dimensions', {ydim_name,ny},...
         'Datatype', 'int32',...
          'FillValue', int32(1e20));
      end
    end

    nccreate(flesper_out, 'time', ...
        'Dimensions', {'time', Inf}, ...
        'Datatype', 'double', ...
        'FillValue', 1e20);

    % Add the units & calendar attributes for time variable:
    dv_ref = datevec(dnmb_ref);
    tref_str = sprintf('days since %04d-%02d-%02d 00:00:00', dv_ref(1), dv_ref(2), dv_ref(3));
    ncwriteatt(flesper_out, 'time', 'units', tref_str);
    ncwriteatt(flesper_out, 'time', 'calendar', 'gregorian');

    % Write data yr by yr
    % format: alk_segment_002(time, nz_segment_002, ny_segment_002, nx_segment_002) 
    itime = 0;
    for iyr = 1:Nyrs_save
      for isegm = 1:4
        segstr = sprintf('_segment_%3.3i',isegm);
        alk_var = sprintf('alk%s',segstr);
        dic_var = sprintf('dic%s',segstr);
        dzalk_var = sprintf('dz_alk%s',segstr);
        dzdic_var = sprintf('dz_dic%s',segstr);
        lat_var   = sprintf('lat%s',segstr);
        lon_var   = sprintf('lon%s',segstr);
        zdim_name = sprintf('nz%s',segstr);
        ydim_name = sprintf('ny%s',segstr);
        xdim_name = sprintf('nx%s',segstr);

        alkT = DATA(iyr).segm(isegm).alk;  % [days, nz, nx]
        dicT = DATA(iyr).segm(isegm).dic; 
        TMM  = DATA(iyr).segm(isegm).TM - dnmb_ref;   % monthly data time stamps, days since ...
        lon_segm = DATA(iyr).segm(isegm).lon;
        lat_segm = DATA(iyr).segm(isegm).lat;
       
        nx  = DATA(iyr).segm(isegm).nx_segm;
        ny  = DATA(iyr).segm(isegm).ny_segm;
        nz  = DATA(iyr).segm(isegm).nz_segm;
        nt  = length(TMM);
      
        dstr1 = datestr(TMM(1)+dnmb_ref, 'yyyy-mm-dd');
        dstr2 = datestr(TMM(end)+dnmb_ref, 'yyyy-mm-dd');
        fprintf('  Segment %i | nx = %i, ny = %i | time: %s to %s\n', isegm, nx, ny, ds1, ds2);

        if nx > 1 && ny > 1
          error(' ERR: one of dim should be singleton: nx=%i, ny=%i\n',nx,ny);
        end

        % Check dimensions:
        assert(isequal(size(alkT), size(dicT)), 'alkT and dicT must have the same size\n');
        [tdim, zdim, hdim] = size(alkT); 
        assert(zdim == nz, 'alkT: zdim=%i should be %i',zdim,nz);
        assert(hdim == max([nx,ny]), 'alkT: hdim=%i should be %i',hdim,max([nx,ny]));
        assert(tdim == nt, 'alkT: tdim=%i should be %i',tdim, nt); 

        % Add a singleton dim: 
        alkT = reshape(alkT, nt, nz, ny, nx); 
        dicT = reshape(dicT, nt, nz, ny, nx); 

        DZ = permute(DZ, [3,2,1]);
        DZ = DZ(1:nt,:,:);
        DZ = reshape(DZ, nt, nz, ny, nx);

        ncwrite(flesper_out, alk_var,   alkT, [itime+1, 1, 1, 1]);
        ncwrite(flesper_out, dic_var,   dicT, [itime+1, 1, 1, 1]);
        ncwrite(flesper_out, dzalk_var, DZ, [itime+1, 1, 1, 1]);
        ncwrite(flesper_out, dzdic_var, DZ, [itime+1, 1, 1, 1]);
        if iyr == 1
          ncwrite(flesper_out, lat_var,   lat_segm, 1);
          ncwrite(flesper_out, lon_var,   lon_segm, 1);

          if isegm == 1 || isegm == 3
            ncwrite(flesper_out, xdim_name, 0:nx-1, 1);
          else
            ncwrite(flesper_out, ydim_name, 0:ny-1, 1);
          end
        end
      end
      itime = itime+nt;
    end

    % Rearrange fields:
    for iiy=2:kYR
      DATA(iiy-1).alk = DATA(iiy).alk;
      DATA(iiy-1).dic = DATA(iiy).dic;
      DATA(iiy-1).TM  = DATA(iiy).TM;
      kYR = 1;
    end

  end    % if kYR

end    % years

        

% Check:
fchck=0;
if fchck
  AMX=ESPER_Mixed([1 2 3 4 5 6 7],[0 0 100;0 0 1000;-150 0 100],...
                    [35 0.5 5 10 20 200;35 0.5 5 10 20 200;32 0.5 5 10 20 200],[1 3 2 4 5 6],...
                    'Equations',[1 16 8],'EstDates',[1980;2002;2030]);

  ANN=ESPER_NN([1 2 3 4 5 6 7],[0 0 100;0 0 1000;-150 0 100],...
                    [35 0.5 5 10 20 200;35 0.5 5 10 20 200;32 0.5 5 10 20 200],[1 3 2 4 5 6],...
                    'Equations',[8],'EstDates',[1980;2002;2030]);


  op=outp_coords(1:4,:);
  pv=pred_vars(1:4,:);
  TLR=ESPER_LIR(requested_vars, op, pv, predictor_types, 'Equations', [8], 'EstDates', est_dates);
  TNN=ESPER_NN(requested_vars, op, pv, predictor_types, 'Equations', [8], 'EstDates', est_dates);
end
 
