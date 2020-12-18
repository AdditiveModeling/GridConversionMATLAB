clear
clc
%Jacob Pessin
%15 December 2020
%Grid Conversion with Consolidation

tic();
global coordinates elements nel nnodes basename np
global psi LM irow icol nzmax phi
%==========================================================================
%Only edit values within these headers
porosity = 0.55;
powder_thick = 30; %Microns
%File location for the last time step of the simulation
basename = '/Users/Jacob/Desktop/result_single_layer/500x_layer_130/0/';
%Desired ouput filename and location for the csv file
fnameCsvOutput = "/Users/Jacob/Desktop/500x_layer_interp.csv";
%File location for the pvd file from the simulation to read out time steps
fnameTimeRead = '/Users/Jacob/Desktop/result_single_layer/500x_layer.pvd';
%File location basename for reading in temperature values at each time step
basename2 = '/Users/Jacob/Desktop/result_single_layer/500x_layer_';

%Coarse Grid Values
Origins = [200,1,17];
Spacing = [0.5,0.5,0.5];
NumPts = [1,150,106];
%NumPts = [1,15,11];

%Fine Grid Spacing
%InterpSpacing = [0.046,0.046,0.046];
InterpSpacing = [0.46,0.46,0.46];
%Time step interval, must be a whole number and designates time steps
%skiped. For example id Inter = 10, then the time steps recorded will be
%1,11,21... and so on.
Inter = 1;

% number of processors
np = 4;
%==========================================================================
% count number of nodes and cells
nnodesG = zeros(1,np);
ncellsG = zeros(1,np);
offset = zeros(1,np);

%This loop automatically gets the number of elements and nodes for each
%partition and writes them in the array nnodesG, ncellsG
for i=0:np-1
    fname = [basename num2str(i) '.vtu'];
    fileID = fopen(fname,'r');
    
    fgetl(fileID);
    fgetl(fileID);
    tline = fgetl(fileID);
    % get number of points
    nnodesl = sscanf(tline,'<Piece NumberOfPoints="%d"');
    nnodesG(i+1) = nnodesl;
    
    % get number of cells
    str2 =['<Piece NumberOfPoints="' num2str(nnodesl) '"'];
    
    str3 = [str2 ' NumberOfCells="%d"'];
    ncellsl = sscanf(tline,str3);
    ncellsG(i+1) = ncellsl;
    fclose(fileID);
end

% This is the total no. of nodes and elements in the entire mesh
nnodes = sum(nnodesG);
nel = sum(ncellsG);

for i=2:np
    offset(i) = offset(i-1) + nnodesG(i-1);
end

% read the coordinates
coordinates = zeros(nnodes,3);
elements = zeros(nel,4);
psi = zeros(nel,1);
phi = zeros(nel,1);

str = '<DataArray type="Float64" Name="coordinates" NumberOfComponents="3" format="ascii">';
count = 0;
fprintf("got here\n");
for i=0:np-1
    fname = [basename num2str(i) '.vtu'];
    fileID = fopen(fname,'r');
    tline = fgetl(fileID);
    while ischar(tline)
        if strcmp(tline,str)
            break;
        end
        tline = fgetl(fileID);
    end
    for k=1:nnodesG(i+1)
        count = count + 1;
        tline = fgetl(fileID);
        x = sscanf(tline,'%f');
        coordinates(count,:) = x(:);
    end
    fclose(fileID);
end
fprintf("read 1\n");

count = 0;
str = '<DataArray type="Int32" Name="connectivity" format="ascii">';
for i=0:np-1
    fname = [basename num2str(i) '.vtu'];
    fileID = fopen(fname,'r');
    tline = fgetl(fileID);
    while ischar(tline)
        if strcmp(tline,str)
            break;
        end
        tline = fgetl(fileID);
    end
    for k=1:ncellsG(i+1)
        count = count + 1;
        tline = fgetl(fileID);
        x = sscanf(tline,'%d');
        elements(count,:) = x(:) + offset(i+1) + 1;
    end
    fclose(fileID);
end
fprintf("read 2\n");


count = 0;
str = '<DataArray type="Float64" Name="Psi1_1" NumberOfComponents="1" format="ascii">';
for i=0:np-1
    fname = [basename num2str(i) '.vtu'];
    fileID = fopen(fname,'r');
    tline = fgetl(fileID);
    while ischar(tline)
        if strcmp(tline,str)
            disp('found');
            break;
        end
        
        tline = fgetl(fileID);
    end
    for k=1:ncellsG(i+1)
        count = count + 1;
        tline = fgetl(fileID);
        psi(count) = str2double(tline);
    end
    fclose(fileID);
end
fprintf("read 3\n");

count = 0;
str = '<DataArray type="Float64" Name="Phi1_1" NumberOfComponents="1" format="ascii">';
for i=0:np-1
    fname = [basename num2str(i) '.vtu'];
    fileID = fopen(fname,'r');
    tline = fgetl(fileID);
    while ischar(tline)
        if strcmp(tline,str)
            break;
        end
        tline = fgetl(fileID);
    end
    for k=1:ncellsG(i+1)
        count = count + 1;
        tline = fgetl(fileID);
        phi(count) = str2double(tline);
    end
    fclose(fileID);
end
fprintf('read 4\n')

[C,ia,ic] = unique(coordinates,'rows','stable');
coordinates = C;
nnodes = length(coordinates);

for i = 1:4*nel
    elements(i) = (ic(elements(i)));
end
fprintf("Fixed boundary nodes\n");
% get boundary nodes

Tol = 1;

count = 0;
for i = 1:nnodes
    z = coordinates(i,3);
    if abs ( z - powder_thick ) < Tol
        count = count + 1;
    end
end

fixnodes = zeros(count,1);
count = 0;
for i = 1:nnodes
    z = coordinates(i,3);
    if abs ( z - powder_thick ) < Tol
        count = count + 1;
        fixnodes(count) = i;
    end
end

% ====================
% assembling ID array
% ====================
ID = ones(1,nnodes);
[ndispl,n] = size(fixnodes);
% prescribed displacements
for i=1:ndispl
    nd = fixnodes(i,1);
    ID(nd) = 0;
end

% Fill ID array
count = 0;
for j=1:nnodes
    if ( ID(j) ~= 0 )
        count = count + 1;
        ID(j) = count;
    end
end
fprintf("assembled ID array\n");

% =================
% Generate LM array
% =================
LM = zeros(4,nel);
for k=1:nel
    for j=1:4
        LM(j,k) = ID(elements(k,j));
    end
end

ndof = max(ID);
% displacement vector
d = zeros(nnodes,1);
fprintf("generated LM array\n");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute Sparcity
nzmax = 0;
for elem=1:nel
    for k=1:4
        i_index = LM(k,elem);
        if (i_index > 0)
            for m=1:4
                j_index = LM(m,elem);
                if (j_index > 0)
                    nzmax = nzmax + 1;
                end
            end
        end
    end
end


irow = zeros(1,nzmax);
icol = zeros(1,nzmax);

count = 0;
for elem=1:nel
    for k=1:4
        i_index = LM(k,elem);
        if (i_index > 0)
            for m=1:4
                j_index = LM(m,elem);
                if (j_index > 0)
                    count = count + 1;
                    irow(count) = i_index;
                    icol(count) = j_index;
                end
            end
        end
    end
end
fprintf("Computed Sparcity\n");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assembling stiffness matrix
% ===========================
%tic();
K = zeros(1,nzmax);
F = zeros(ndof,1);
count = 0;
for i=1:nel %loop over elements
    xe = coordinates(elements(i,:),:);
    Psie = psi(i);
    [ke,fe] = weakform(xe,Psie,porosity);
    for j=1:4
        i_index = LM(j,i);
        if (i_index > 0)
            F(i_index) = F(i_index) + fe(j);
            for k=1:4
                j_index = LM(k,i);
                if (j_index > 0)
                    count = count + 1;
                    K(count) = K(count) + ke(j,k);
                end
            end
        end
    end
end
fprintf("assembled stiffness matrix\n");
% ensamble sparse matrix
M = sparse(irow,icol,K,ndof,ndof);
fprintf("ensambled sparse matrix\n")
% solve system of equations
a_bar = M\F;
fprintf("Solved system of equations\n");
% change small values (1.0e-12) to zero
[a_bar] = changeZero(a_bar);
fprintf("Change small values to 0\n");
% corrector phase
for i=1:nnodes
    index = ID(i);
    if ( index ~= 0 )
        d(i) = a_bar(index);
        if d(i) > powder_thick
            d(i) = porosity*powder_thick;
            %d(i) = 0;
        end
        if isnan(d(i)) == 1
            d(i) = 0;
            %d(i) = porosity*powder_thick;
        end
    end
end
fprintf("finished corrector phase\n");
coordinates(:,3) = coordinates(:,3) + d;

%Calculate coarse grid to be interpolated onto
xmax = (NumPts(1)*Spacing(1)+Origins(1)-Spacing(1));
ymax = (NumPts(2)*Spacing(2)+Origins(2)-Spacing(2));
zmax = (NumPts(3)*Spacing(3)+Origins(3)-Spacing(3));
xGrid = Origins(1):Spacing(1):xmax;
yGrid = Origins(2):Spacing(2):ymax;
zGrid = Origins(3):Spacing(3):zmax;
[x,y,z] = meshgrid(xGrid,yGrid,zGrid);

%Reduce Element list to search
count = 1;
for i = 1:length(elements)
    if sum((coordinates(elements(i,:),1) - xmax) < 4) == 4 &&...
            sum((coordinates(elements(i,:),2) - ymax) < 4) == 4 &&...
            sum((coordinates(elements(i,:),3) - zmax) < 4) == 4 &&...
            sum((coordinates(elements(i,:),1) - Origins(1)) > -4) == 4 &&...
            sum((coordinates(elements(i,:),1) - Origins(2)) > -4) == 4 &&...
            sum((coordinates(elements(i,:),1) - Origins(3)) > -4) == 4
        elemReduced(count,:) = elements(i,:);
        count = count + 1;
    end
end
disp('Finished reducing element search list')

%Find Elements associated with each grid point
elemAssoc = zeros(length(xGrid)*length(yGrid)*length(zGrid),4);
for l = 1:length(xGrid)
    for i = 1:length(yGrid)
        for j = 1:length(zGrid)
            for k = 1:length(elemReduced)
                if inhull([xGrid(l),yGrid(i),zGrid(j)],coordinates(elemReduced(k,:),:)) == 1
                    elemAssoc(l+(i-1)*length(xGrid)+(j-1)*length(yGrid)*length(xGrid),:) = elemReduced(k,:);
                    break
                end
            end
            %disp(j)
        end
        disp(i)
    end
end
disp('Finished finding which element each grid point is in')

%Loop through timesteps to interpolate temperature data onto grid
newTemp = zeros(130/Inter,length(elemAssoc));
count2 = 1;
for j = 1:Inter:130
    basename3 = [basename2 num2str(j) '/0/'];
    temp = zeros(nnodes,1);
    str = '<DataArray type="Float64" Name="temp_old" NumberOfComponents="1" format="ascii">';
    count = 0;
    fprintf("got here\n");
    for i=0:np-1
        fname = [basename3 num2str(i) '.vtu'];
        fileID = fopen(fname,'r');
        tline = fgetl(fileID);
        while ischar(tline)
            if strcmp(tline,str)
                break;
            end
            tline = fgetl(fileID);
        end
        for k=1:nnodesG(i+1)
            count = count + 1;
            tline = fgetl(fileID);
            t = sscanf(tline,'%f');
            temp(count,:) = t(:);
        end
        fclose(fileID);
    end
    fprintf("read temp\n");
    for i = 1:length(xGrid)
        for m = 1:length(yGrid)
            for k = 1:length(zGrid)
                l = i+(m-1)*length(xGrid)+(k-1)*length(yGrid)*length(xGrid);
                x1 = [coordinates(elemAssoc(l,1),1),coordinates(elemAssoc(l,2),1),...
                    coordinates(elemAssoc(l,3),1),coordinates(elemAssoc(l,4),1)];
                x2 = [coordinates(elemAssoc(l,1),2),coordinates(elemAssoc(l,2),2),...
                    coordinates(elemAssoc(l,3),2),coordinates(elemAssoc(l,4),2)];
                x3 = [coordinates(elemAssoc(l,1),3),coordinates(elemAssoc(l,2),3),...
                    coordinates(elemAssoc(l,3),3),coordinates(elemAssoc(l,4),3)];
                %         xx = x'; xq1 = xx(l);
                %         yy = y'; xq2 = yy(l);
                %         zz = z'; xq3 = zz(l);
                xq1 = xGrid(i);
                xq2 = yGrid(m);
                xq3 = zGrid(k);
                V = [temp(elemAssoc(l,1)),temp(elemAssoc(l,2)),...
                    temp(elemAssoc(l,3)),temp(elemAssoc(l,4))];
                newTemp(count2,l) = griddata(x1,x2,x3,V,xq1,xq2,xq3);
            end
        end
    end
    count2 = count2 + 1;
end
disp('Interpolated the temperature onto the coarse grid')

fileID = fopen(fnameTimeRead,'r');
fgetl(fileID);
fgetl(fileID);
fgetl(fileID);
for j = 1:130
    tline = fgetl(fileID);
    time(j) = sscanf(tline,'      <DataSet timestep="%f"');
end
fclose(fileID);
timeRed = time(1:Inter:130);
disp('Finished reading timesteps');

xqGrid = Origins(1):InterpSpacing(1):(NumPts(1)*Spacing(1)+Origins(1)-Spacing(1));
yqGrid = Origins(2):InterpSpacing(2):(NumPts(2)*Spacing(2)+Origins(2)-Spacing(2));
zqGrid = Origins(3):InterpSpacing(3):(NumPts(3)*Spacing(3)+Origins(3)-Spacing(3));
[xq,yq,zq] = meshgrid(xqGrid,yqGrid,zqGrid);

InterpNumPts = [length(xqGrid),length(yqGrid),length(zqGrid)];

fid = fopen (fnameCsvOutput, 'w');
fprintf(fid,'Regular Grid Temperature History for Simulation 500x_layer_interp\n');
fprintf(fid,'%s,%s,%s\n','X Origin','Y Origin','Z Origin');
fprintf(fid,'%d,',Origins); fprintf(fid,'\n');
fprintf(fid,'%s,%s,%s\n','X Spacing','Y Spacing','Z Spacing');
fprintf(fid,'%d,',InterpSpacing); fprintf(fid,'\n');
fprintf(fid,'%s,%s,%s\n','Num Pts X','Num Pts Y','Num Pts Z');
fprintf(fid,'%d,',InterpNumPts); fprintf(fid,'\n');
fprintf(fid,'%s,%s,%s,%s\n','Time (seconds)','Temperature at Regular Grid Points (nested x','then y','then z list)');

count = 1;
for i = 10:(length(timeRed)+9)    
    if NumPts(1) > 1 && NumPts(2) > 1 && NumPts(3) > 1
        TempShape = reshape(newTemp(count,:),[NumPts(2),NumPts(1),NumPts(3)]);
        Tempq = interp3(x,y,z,TempShape,xq,yq,zq);
        [l,r,c] = size(Tempq);
        TempqLin = reshape(Tempq,[1,l*r*c]);
    else
        TempShape = reshape(newTemp(count,:),[NumPts(3),NumPts(2)]);
        [y2,z2] = meshgrid(yGrid,zGrid);
        [yq2,zq2] = meshgrid(yqGrid,zqGrid);
        Tempq = interp2(y2,z2,TempShape,yq2,zq2);
        [r,c] = size(Tempq);
        TempqLin = reshape(Tempq,[1,r*c]);
    end
    
    fprintf(fid,'%d,%d,',timeRed(count),TempqLin); fprintf(fid,'\n');
    count = count + 1;
end
fclose(fid);

elapsed_time = toc();
fprintf('Elapsed time (s) to postprocess data = %i\n', elapsed_time);