function Vm=extractFrame(filename_bin,frame,filename_geom)
%Vm=extractFrame(filename_bin,frame) read the binary file with name filename_bin and
%extract the indicated frame. plotFrame search for a mat
%file in the same path and with the same name containing simulation
%metadata (geometry adn fibrosis). Otherwise the mat file filename can be
%passed as third argument to the extractFrame function.

fid=fopen(filename_bin,'r');

if nargin<3
    filename_geom=filename_bin;
end


load(filename_geom,"ind_in");
sizeVm=nnz(ind_in);
fseek(fid,0,'eof');
Nframe=ftell(fid)/4/sizeVm;

if frame>Nframe || frame<1
    fclose(fid);
    error('Selcted frame does not exists');
end

fseek(fid,4*sizeVm*(frame-1),'bof');
Vm=fread(fid,sizeVm,'single'); 

end


