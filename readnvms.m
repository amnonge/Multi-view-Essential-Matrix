
function [groundKs,groundRs,groundTs,files]=readnvms(file)

fid = fopen(sprintf('%s',file));
tline = fgets(fid);
tline = fgets(fid);
tline = fgets(fid);
numOfCameras=str2num(tline);
tline = fgets(fid);

cams=zeros(numOfCameras,8);
files=cell(numOfCameras,1);
for i=1:numOfCameras
    C = strsplit(tline);
%     cc=C(1);
    fileName=C{1};
    files{i}=fileName;
    cams(i,:)=[str2double(C{2}) str2double(C{3}) str2double(C{4}) str2double(C{5}) str2double(C{6}) str2double(C{7}) str2double(C{8}) str2double(C{9})];
    tline = fgets(fid);
end

fclose(fid);

groundKs=cell(numOfCameras,1);
groundRs=cell(numOfCameras,1);
groundTs=cell(numOfCameras,1);

for i=1:length(files)
    curFile=files{i};
    cam=cams(i,2:8);
    focal=cams(i,1);
%     im=imread(curFile);
%     groundKs{i}=[focal 0 30size(im,2)/2;0 focal size(im,1)/2;0 0 1];
    curQ=cam([1 2 3 4]);
    curt=cam(5:7);
    groundRs{i}=quat2rotm(curQ)';
    groundTs{i}=curt';
    

end

end