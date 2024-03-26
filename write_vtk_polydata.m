function write_vtk_polydata(filename,node,elem,scalar,colorscale)
% write_vtk_polydata(filename,node,elem,scalar,colorscale) write a surface 
% mesh defined by node, elem and scalar values to a vtk polydata file with
% name filename. A colorscale should be also defined for color representation of the
% scalar value.
fid=fopen(filename,'w+');
fprintf(fid, '# vtk DataFile Version 4.1\n');
fprintf(fid, 'PatientData from Matlab\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET POLYDATA\n');
fprintf(fid, ['POINTS ' num2str(length(node)) ' float\n']);
fprintf(fid,[repmat('%0.5f ', 1,9) '\n'],node');
fprintf(fid,'\n');
fprintf(fid,['POLYGONS ' num2str(length(elem)) ' ' num2str(length(elem)*4) '\n']);
fprintf(fid,'3 %u %u %u \n',elem'-1);
fprintf(fid,['POINT_DATA ' num2str(length(node)) '\n' ]);
fprintf(fid,'SCALARS act_times float \n');
fprintf(fid,'LOOKUP_TABLE lookup_table \n');
fprintf(fid,[repmat('%0.5f ', 1,9) '\n'],scalar);
fprintf(fid,'\n');
fprintf(fid,['LOOKUP_TABLE lookup_table ' num2str(size(colorscale,1)) '\n']);
fprintf(fid,'%0.5f %0.5f %0.5f %0.5f \n',[colorscale ones(size(colorscale,1),1)]');
fclose(fid);
end
