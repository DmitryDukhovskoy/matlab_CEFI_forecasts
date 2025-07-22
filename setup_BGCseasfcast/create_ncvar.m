function create_ncvar(ncfile, varname, dims, dtype)
  % delete(flesper_out);
  %
  % Default for optional flag
  if nargin < 4
    dtype = 'double';
  end

  % Check if variable already exists
  if isfile(ncfile)
    vars = {ncinfo(ncfile).Variables.Name};
    if any(strcmp(varname, vars))
      fprintf('%s already exists\n', varname);
      return;
    end
  end

  % Create variable
  if strcmp(dtype, 'int32')
    nccreate(ncfile, varname, ...
       'Dimensions', dims, ...
       'Datatype', 'int32', ...
       'Format', 'netcdf4');
  else
    nccreate(ncfile, varname, ...
       'Dimensions', dims, ...
       'Datatype', 'double', ...
       'FillValue', 1e20, ...
       'Format', 'netcdf4');
  end

  fprintf('     %s created\n', varname);
end

