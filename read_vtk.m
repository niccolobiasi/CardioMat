function VTKData=read_vtk(filename)
% read_vtk read vtk polydata file format of ADAS shells
% VTKData=read_vtk(filename) produces a structure VTKData by reading the 
% filename vtk file with the following numeric fields:
% -points (NX1): contains the coordinates of the N nodes in the vtk file
% -elements (NelX3): contains the nodes indices composing each of the Nel
% elements
% -scalar(NX1): contains the N scalar values for each node
%
% NOTE: this function is only inteded to read ADAS vtk shells and not any
% other vtk file that may have a different structures. Adjustments are
% needed for other vtk files.
%
VTKData.header = cell(4,1);

fid = fopen(filename,'r');

if fid==-1
    error("file not found")
end

% 4 lines of header
VTKData.header{1} = fgetl(fid);
if (strcmp(VTKData.header{1}(3:5),'vtk') ~= 1)
    disp('Verify vtk file type and header structure');
    return;
end
VTKData.header{2} = fgetl(fid);
VTKData.header{3} = fgetl(fid);
VTKData.header{4} = fgetl(fid);
VTKData.gridType = VTKData.header{4}(9:end);
% flags for searching a title line to find where the numbers are
flag1 = 0;
flag2 = 0;
str='';

while isempty(str)
    str = fgetl(fid);
end

if (strcmp(str(1:6),'POINTS') == 1)
    % find the end of the number
    for i = 7:length(str)
        if (str(i) ~= ' ' && flag1 == 0)
            flag1 = 1;
        end
        if (str(i) == ' ' && flag1 == 1)
            numPoints = str2double(str(7:i-1));
            break;
        end
    end
    % initialize the data structure
    VTKData.points = zeros(numPoints,3);
    kk=1;
    while kk<=numPoints
        str = fgetl(fid);
        data_line=sscanf(str,'%f');
        num_data=length(data_line);
        num_points_line=num_data/3;
        data_array= (reshape(data_line,[3,num_points_line]))';
        VTKData.points(kk:kk+num_points_line-1,:) = data_array;
        kk=kk+num_points_line;
    end
else
    error('Point field not found')
end

str='';
while isempty(str)
    str = fgetl(fid);
end
if (strcmp(str(1:8),'POLYGONS') == 1)
    for i = 9:length(str)
        if (str(i) ~= ' ' && flag2 == 0)
            flag2 = 1;
        end
        if (str(i) == ' ' && flag2 == 1)
            numCells = str2double(str(9:i-1));
            break;
        end
    end


    VTKData.elements=zeros(numCells,3);
    for kk=1:numCells
        str = fgetl(fid);
        data_line=sscanf(str,'%f');
        points_cell=data_line(2:end)+1;
        VTKData.elements(kk,:)=points_cell;
    end

else
    error('Polygons field not found')
end

str='';
while isempty(str)
    str = fgetl(fid);
end

if strcmp(str(1:10),'POINT_DATA')
    numPoints=sscanf(str(11:end),'%f');
    str=fgetl(fid);
    while isempty(sscanf(str,'%f'))
        str=fgetl(fid);
    end
    VTKData.scalar = zeros(numPoints,1);
    kk=1;
    while kk<=numPoints
        data_line=sscanf(str,'%f');
        num_data=length(data_line);
        VTKData.scalar(kk:kk+num_data-1,:) = data_line';
        kk=kk+num_data;
        str = fgetl(fid);
    end
else
    warning('No data points found')
end

fclose(fid);

end
