function extracted_lat= extractLat(lat,type)
% extracted_lat= extractLat(lat,type) extract a specific lat from the cell
% array outputted by the computeLats function. If the type input is an
% integer extract_lat returns a numeric array containing for each voxel
% the activation identified by type. Type can also be 'first' or 'last'
% string to output the first or the last activation respectively.

    fun=@(x) single_fun(x,type);
    extracted_lat=cellfun(fun,lat);
end

function lat_out=single_fun(lat,type)
if isempty(lat)
    lat_out=NaN;
elseif isnumeric(type)
    lat_out=lat(type);
elseif ischar(type)
    if strcmp(type,'first')
        lat_out=lat(1);
    elseif strcmp(type,'last')
        lat_out=lat(end);
    end
end
end
