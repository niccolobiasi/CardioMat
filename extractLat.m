function extracted_lat= extractLat(lat,type,start_t,end_t)
% extracted_lat= extractLat(lat,type) extract a specific lat from the cell
% array outputted by the computeLats function. If the type input is an
% integer extract_lat returns a numeric array containing for each voxel
% the activation identified by type. Type can also be 'first' or 'last'
% string to output the first or the last activation respectively.

if nargin<3
    start_t=0;
end

if nargin<4
    end_t=inf;
end


fun=@(x) single_fun(x,type,start_t,end_t);
extracted_lat=cellfun(fun,lat);
end

function lat_out=single_fun(lat,type,start_t,end_t)
if isempty(lat)
    lat_out=nan;
elseif isnumeric(type)
    if length(lat)>=type
        lat_out=lat(type);
    else
        lat_out=nan;
    end
elseif ischar(type)
    if strcmp(type,'first')
        tmp=lat(lat<end_t & lat>start_t);
        if isempty(tmp)
            lat_out=nan;
        else
            lat_out=tmp(1);
        end
    elseif strcmp(type,'last')
        tmp=lat(lat<end_t & lat>start_t);
        if isempty(tmp)
            lat_out=nan;
        else
            lat_out=tmp(end);
        end
    end
end
end
