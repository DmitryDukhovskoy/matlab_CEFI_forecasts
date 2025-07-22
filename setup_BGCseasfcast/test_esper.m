% The code has not been fully tested 
%
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

nz=75;
nt=13;
nx=120;
ny=1;
t2D = rand(nt,nz);
t4D = reshape(t2D, [nt, nz, 1, 1]);
t4D = permute(t4D, [4 3 2 1]); 
lat_segm = rand(nx,1);
lat_var = 'x_segm';

flnm_test = sprintf('test_esper_SPEARmnth_%i-%i.nc',yr_start,yr_end);
fltest_out = fullfile(pthesper, flnm_test);

if isfile(fltest_out)
  delete(fltest_out);
end
fprintf('Saving fields --> %s\n',fltest_out);

% Note Matlab saves to betcdf in reverse order 
nccreate(fltest_out, 'T4D', ...
   'Dimensions', {'xdim',nx, 'ydim',ny, 'zdim',nz, 'time',Inf}, ...
   'Datatype','double',...
   'FillValue', 1e20,...
   'Format', 'netcdf4');

nccreate(fltest_out, lat_var, ...
 'Dimensions', {'xdim',nx},...
 'Datatype', 'double',...
 'FillValue', 1e20,...
 'Format', 'netcdf4');


nccreate(fltest_out, 'time', ...
    'Dimensions', {'time', Inf}, ...
    'Datatype', 'double', ...
    'FillValue', 1e20, ...
    'Format', 'netcdf4');

irec = 0;
ncwrite(fltest_out, 'T4D', t4D, [1 1 1 irec+1]);
if irec == 0
  ncwrite(fltest_out, lat_var, lat_segm);
end

irec = irec + nt;
ncwrite(fltest_out, 'T4D', t4D, [1 1 1 irec+1]);

        

