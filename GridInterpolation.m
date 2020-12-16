clear
clc
%Jacob Pessin
%Grid Conversion Interpolation
%InterpSpacing = [0.046,0.046,0.046];
InterpSpacing = [0.46,0.46,0.46];

filename = "/Users/Jacob/Desktop/500x_layer.csv";
Origins = csvread(filename,2,0,[2,0,2,2]);
Spacing = csvread(filename,4,0,[4,0,4,2]);
NumPts = csvread(filename,6,0,[6,0,6,2]);
yGrid = Origins(2):Spacing(2):(NumPts(2)*Spacing(2)+Origins(2)-Spacing(2));
zGrid = Origins(3):Spacing(3):(NumPts(3)*Spacing(3)+Origins(3)-Spacing(3));
[y,z] = meshgrid(yGrid,zGrid);

fname = "/Users/Jacob/Desktop/500x_layer_interp_test.csv";
fid = fopen (fname, 'w');
fprintf(fid,'Regular Grid Temperature History for Simulation 500x_layer_interp\n');
fprintf(fid,'%s,%s,%s\n','X Origin','Y Origin','Z Origin');
fprintf(fid,'%d,',Origins); fprintf(fid,'\n');
fprintf(fid,'%s,%s,%s\n','X Spacing','Y Spacing','Z Spacing');
fprintf(fid,'%d,',InterpSpacing); fprintf(fid,'\n');
fprintf(fid,'%s,%s,%s\n','Num Pts X','Num Pts Y','Num Pts Z');
fprintf(fid,'%d,',NumPts); fprintf(fid,'\n');
fprintf(fid,'%s,%s,%s,%s\n','Time (seconds)','Temperature at Regular Grid Points (nested x','then y','then z list)');

for i = 10:137
    Time = csvread(filename,i,0,[i,0,i,0]);
    Temp = csvread(filename,i,1,[i,1,i,NumPts(2)*NumPts(3)]);
    Temp = reshape(Temp,[NumPts(3),NumPts(2)]);
        
    yqGrid = Origins(2):InterpSpacing(2):(NumPts(2)*Spacing(2)+Origins(2)-Spacing(2));
    zqGrid = Origins(3):InterpSpacing(3):(NumPts(3)*Spacing(3)+Origins(3)-Spacing(3));
    [yq,zq] = meshgrid(yqGrid,zqGrid);
    Tempq = interp2(y,z,Temp,yq,zq);
    [r,c] = size(Tempq);
    TempqLin = reshape(Tempq,[1,r*c]);
    
    fprintf(fid,'%d,%d,',Time,TempqLin); fprintf(fid,'\n');
end
fclose(fid);
figure
subplot(1,2,1)
surf(Temp,y,z)
subplot(1,2,2)
surf(Tempq,yq,zq)