 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is associated with the paper titled "Cytoskeletal activation of
% NHE1 regulates mechanosensitive cell volume adaptation and proliferation"
% by Ni, et al. 
% The code is can be used for academic and research purposes only.
% For inquiries, please contact Qin Ni or Sean X. Sun at Johns Hopkins
% University. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constants
% AimFolder is the folder where all the .tiff files are saved into
AimFolder = pwd;
% Timelengthforgoodcells = (60+120)/10+1; % The length of cells viedos to be considered as good
cd(AimFolder) % Go to the AimFolder Automatically
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumberofTiles = size(dir,1)-2;
mkdir(strcat(AimFolder,'\','Combine'))

% Parameters for cell cropping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chHeight = 12; % This the estimated channel height only for cell cropping
Filter0 =1.7; % This is used to locate the cells
hch = 280; % half the length of the largest crop cell area
Dilate = 70;
mnp = 1500; % During mask the minimum number of connected pixels to be considered as a cell


%% Run for different folders
for tile = 1:NumberofTiles
    CurrentTile = tile    
    Foldername = strcat(AimFolder,'\',string(tile));
    cd(Foldername)
    rector = dir('*.tif');
    numCh = 1; % 
    R = length(rector)/numCh;
    DataOri = cell(R,numCh+2); % Original Images
    DataPos = cell(R,1);
    Location = cell(1,R);
    CroppingEdge = {}; % Used to crop cells
    OrderofCrop = {};
    Distance = {}; % Clear before next tile
    DistanceBackward = {}; % Celar before next tile
    % Extracting single cells Automatically
    f = 1;
    SnapShotf = cell(1,numCh);
    for i = 1:R
        for j =  1: numCh
            SnapShotf{j} = imread(rector((i-1)*numCh+j).name);
        end
        Vchannel = double(SnapShotf{1});
        sizei = size(Vchannel);
        [c,r] = meshgrid(linspace(1,sizei(2),sizei(2)),linspace(1,sizei(1),sizei(1))); 
        %VchannelFit = Vchannel(700:1500,500:2200);
        VchannelFit = Vchannel(700:1500,500:1500);
        sizeIm = size(VchannelFit);
        [colN,rowN] = meshgrid(linspace(1,sizeIm(2),sizeIm(2)),linspace(1,sizeIm(1),sizeIm(1)));
        % box = logical(gt(colN,max(colN(:))-10)+lt(colN,11)+gt(rowN,max(rowN(:))-10)+lt(rowN,11));
        box = logical(gt(colN,max(colN(:))-20)+lt(colN,21)+gt(rowN,max(rowN(:))-20)+lt(rowN,21)); % Thicker Edge
%         planeVfit = fit([colN(:),rowN(:)],VchannelFit(:),'poly22','weights',box(:));
%         planeV = planeVfit.p00 + planeVfit.p10*c + planeVfit.p01*r + planeVfit.p20*c.^2 + planeVfit.p11*r.*c+planeVfit.p02*r.^2;
        planeVfit = fit([colN(:),rowN(:)],VchannelFit(:),'poly11','weights',box(:));
        planeV = planeVfit.p00 + planeVfit.p10*c + planeVfit.p01*r;        
        normV = chHeight*(1-Vchannel./planeV); %%Since for channel, cell is the dark part!!!!!
        n = double(normV > Filter0);
        Mask1 = bwareaopen(n,mnp);
        L = bwlabel(Mask1,8);
%         figure(1)
%         subplot(1,3,1)
%         imagesc(Vchannel)
%         subplot(1,3,2)
%         imagesc(normV)         
%         subplot(1,3,3)
%         imagesc(L)
%         drawnow
        CC = bwconncomp(Mask1);
        for k = 1:size(CC.PixelIdxList,2)
            [Y,X] = ind2sub(size(L),find(L==k)); % So X is horizontal(X's column), Y is vertical(Y's row)
            if max(X)-min(X)<2*(hch-Dilate) & max(Y)-min(Y)<2*(hch-Dilate)
                Xave = mean(X);
                Yave = mean(Y);
                if min(Y)-Dilate>1 & max(Y)+Dilate<sizei(1)& min(X)-Dilate>1 & max(X)+Dilate<sizei(2)...
                        & sum(Mask1(max(Y)+Dilate,min(X)-Dilate:max(X)+Dilate)) == 0 ... %not make the another cell get into the box. Top line all points equal to zero
                        & sum(Mask1(min(Y)-Dilate,min(X)-Dilate:max(X)+Dilate)) == 0 ... % Bottom line
                        & sum(Mask1(min(Y)-Dilate:max(Y)+Dilate,min(X)-Dilate)) == 0 ...
                        & sum(Mask1(min(Y)-Dilate:max(Y)+Dilate,max(X)+Dilate)) == 0 ...            
                    Location{k,i} = [Xave Yave];
                    CroppingEdge{k,i} = [min(X)-Dilate max(X)+Dilate min(Y)-Dilate max(Y)+Dilate];
                end
            end
        end
%        CurrentTime = i
    end
    Location1 = Location;
    % Make empty arrays in Location be [0 0]
    for i = 1:size(Location,1)
        for j = 1:size(Location,2)
            a = isempty(Location{i,j});
            if a == 1
                Location{i,j} = [0 0];
                CroppingEdge{i,j} = [0 0 0 0];
            end
        end
    end
    % Distance Forward
    for c = 1:size(Location,1)
        for t = 1:R-1
            for c2 = 1:size(Location,1)
                Distance{c,t}(c2,1) = (Location{c,t}(1,1)-Location{c2,t+1}(1,1))^2+(Location{c,t}(1,2)-Location{c2,t+1}(1,2))^2;
                %Distance{c,t}(c2,1) = (Location{c,t}(1,1)-Location{c2,t+1}(1,1))^2;
            end
        end
    end
    % Distance Backward
    for c = 1:size(Location,1)
        t = size(Location,2)-1;
        while t >= 1
            for c2 = 1:size(Location,1)
                DistanceBackward{c,t+1}(c2,1) = (Location{c,t+1}(1,1)-Location{c2,t}(1,1))^2+(Location{c,t+1}(1,2)-Location{c2,t}(1,2))^2;
            end
            t = t - 1;
        end
    end
    %
    Forward = {};
    ForwardCrop = {};
    Forward(:,1) = Location(:,1);
    ForwardCrop(:,1) = CroppingEdge(:,1);
    c1 = 1;
    while c1 <= size(Location,1)
        c = c1;
        for t = 1:R-1
            minimumi = min([Distance{c,t}]);
            if minimumi ~= 0 
                temp = Distance{c,t};
                [rowi,columni]=find(temp == minimumi(1));
                Forward{c1,t+1} = Location{rowi,t+1};
                ForwardCrop{c1,t+1} = CroppingEdge{rowi,t+1};
            else
                Forward{c1,t+1} = [0 0];
                ForwardCrop{c1,t+1} = [0 0 0 0];
                rowi = 1;
            end
            c = rowi;
        end
        c1 = c1 + 1;
    end
    %
    Backward = {};
    Backward(:,size(Location,2)) = Location(:,size(Location,2));
    c1 = 1;
    while c1 <= size(Location,1)
        t = size(Location,2);
        c = c1;
        while t > 1
            minimum = min([DistanceBackward{c,t}]);
            if minimum ~= 0 
                temp = DistanceBackward{c,t};
                [rowi,columni]=find(temp== minimum(1));
                Backward{c1,t-1} = Location{rowi,t-1};
            else
                Backward{c1,t-1} = [0 0];
                rowi = 1;  
            end
            t = t - 1;
            c = rowi;
        end
        c1 = c1+1;
    end
    % Move Cell trials with positions of [0 0] from Forward and Backward
    i = size(Forward,1);
    while i >= 1
        iffindzero = 0;
        for j = 1:size(Forward,2)
            a = Forward{i,j}(1,1);
            if a == 0
                iffindzero = 1;
                break
            end
        end

        if iffindzero == 1
            Forward(i,:) = [];
            ForwardCrop(i,:) = [];
        end
        i = i-1;
    end

    i = size(Backward,1);
    while i >= 1
        iffindzero = 0;
        for j = 1:size(Backward,2)
            a = Backward{i,j}(1,1);
            if a == 0
                iffindzero = 1;
                break
            end
        end

        if iffindzero == 1
            Backward(i,:) = [];
        end
        i = i-1;
    end
    % c is the number
    CropPositions = {};
    FinalCropEdge = {};
    iffindsame = 0;
    Timelengthforgoodcells = size(Forward,2);
    for i = 1: size(Forward,1)
        for j = 1: size(Backward,1)
            a = Forward(i,1:Timelengthforgoodcells); % Time point 1-20
            b = Backward(j,1:Timelengthforgoodcells); % Time
            if isequal(a,b)
                iffindsame(i,1) = 1;
                break
            end
        end
    end
    j = 1;
    for i = 1:size(iffindsame,1)
        a = iffindsame(i,1);
        if a == 1
            CropPositions(j,:) = Forward(i,:);
            FinalCropEdge(j,:) = ForwardCrop(i,:);
            j = j + 1;
        end
    end
    % Crop Cells According to the CropPositions
    SnapShotf = cell(1,numCh);
    pic = {};
    for i = 1:R
        for j =  1: numCh
            SnapShotf{i,j} = imread(rector((i-1)*numCh+j).name);
        end 
    end
    for i = 1: size(CropPositions,1) % Available Cells
        for j = 1: size(CropPositions,2) % Time Points
            for k = 1:numCh
                pic = SnapShotf{j,k};
                cellthis = pic(FinalCropEdge{i,j}(1,3):FinalCropEdge{i,j}(1,4),FinalCropEdge{i,j}(1,1):FinalCropEdge{i,j}(1,2));
                DataOri{j,4*(i-1)+k} = double(cellthis);
            end
            DataOri{j,4*(i-1)+3} = CropPositions{i,j}; % This is the mean maybe "center of mass" of the cell
            DataOri{j,4*(i-1)+4} = FinalCropEdge{i,j};
        end
    end
    dataOutput = DataOri;
    savename = strcat(AimFolder,'\Combine','\',string(tile),'rawSingleCells.mat');
    save(savename,'dataOutput','-v7.3')
    NumberofCellsThisTile = size(DataOri,2)/4

end