 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is associated with the paper titled "Cytoskeletal activation of
% NHE1 regulates mechanosensitive cell volume adaptation and proliferation"
% by Ni, et al. 
% The code is can be used for academic and research purposes only.
% For inquiries, please contact Qin Ni or Sean X. Sun at Johns Hopkins
% University. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% combine files from different tiles
clc; clear;
rector = dir('*.mat');
load(rector(1).name);
AlldataOutput = dataOutput;
for i = 2:size(rector,1)
    a = size(AlldataOutput,2);
    load(rector(i).name)
    b = size(dataOutput,2);
    AlldataOutput(:,a+1:a+b) = dataOutput;
end
% Get rid of empty cells
i = size(AlldataOutput,2);
while i >= 1
    a = isempty(AlldataOutput{1,i});
    if a == 1
        AlldataOutput(:,i) = [];
    end
    i = i-1;
end
% select
f = 1;
imax = size(AlldataOutput,2)/3;
i = 1;
while i <= imax
    Data = AlldataOutput(1:end,3*(i-1)+1);
    set(figure,'position',[100 100 1600 1000])
    for j = 1:size(Data,1)
        % subplot(4,5,j)
        subplot(4,9,j)
        imagesc(Data{j,1})
    end
    w = waitforbuttonpress;
    if w == 1
    else
        Cells(1,f) = i;
        f = f + 1;
    end
    i = i + 1;
    close all
end
% Move not good cells
if Cells ~=0
    CellNumber = sort(Cells,'ascend'); % Tiles seperated with "Space"
    i = size(CellNumber,2);
    while i >= 1
        a = CellNumber(1,i);
        AlldataOutput(:,3*a-2:3*a) = [];
        i = i - 1;
    end
end
% Change the format and Save
dataOutput = {};
for i = 1: size(AlldataOutput,2)/3
    for j = 1:size(AlldataOutput,1)
        dataOutput{j,i}(:,:,1) = AlldataOutput{j,3*i-2};
        %dataOutput{j,i}(:,:,2) = AlldataOutput{j,3*i-1};
    end
end
save('rawSingleCells.mat','dataOutput','-v7.3')

%% FXm calculations
%clc;clear all;close
%load('rawSingleCells');
DataOri = dataOutput;
chHeight =12;
%
numCh = 1;
Filter1 = 0.3; %The filter used to find the cell mask
Filter2 = 0.4; %The filter used to crop cell out from background
imSizeMax = 512;

% Image Processing 
DataPro255 = cell(size(DataOri,2),1);
planefitstd = nan(size(DataOri));

for celli = 1 : size(DataOri,2)
    celli
    DataPro255{celli} = cell(size(DataOri,1),numCh+2);
    
    for i = 1 : size(DataOri,1)
        
        if isempty(DataOri{i,celli})
            break
        end
        % Get Cell Mask by Polynomial fitting
        %Vchannel = DataOri{i,numCh} ;
        %minus bg computed outside the channel
        Vchannel = DataOri{i,celli}(:,:,numCh) -300;
        
        sizeIm = size(Vchannel);
        [colN,rowN] = meshgrid(linspace(1,sizeIm(2),sizeIm(2)),linspace(1,sizeIm(1),sizeIm(1)));
        box = double(logical(gt(colN,max(colN(:))-10)+lt(colN,11)+gt(rowN,max(rowN(:))-10)+lt(rowN,11)));
        planeVfit = fit([colN(:),rowN(:)],Vchannel(:),'poly22','weights',box(:));
        planeV = planeVfit.p00 + planeVfit.p10*colN+planeVfit.p01*rowN+planeVfit.p20*colN.^2+planeVfit.p11*rowN.*colN+planeVfit.p02*rowN.^2;
        planefitstd(i,celli) = std(planeV(:));
        normV = chHeight*(1-Vchannel./planeV);
        

        n = double(normV > Filter1);
        Mask = bwareaopen(n,500);
        Mask = imclose(Mask,strel('disk',50));
        dilatesize = 10;
        MaskNew = imdilate(Mask,strel('disk',dilatesize));

        segmentedV = normV.*MaskNew;
        % Padding
        canvasV = zeros(imSizeMax,imSizeMax);
        [s1,s2] = size(segmentedV);
        startRow = max(ceil((imSizeMax-1-s1)/2),1);
        startCol = max(ceil((imSizeMax-1-s2)/2),1);
        canvasV(startRow:startRow+s1-1,startCol:startCol+s2-1) = segmentedV;
        canvasV = imgaussfilt(canvasV,6); 
        canvasV(lt(canvasV,0)) = 0;
        % Calculate volume
        volume = sum(canvasV(:))*0.325^2;% pixel size = 0.335
        
        
        DataPro255{celli}{i,numCh+1} = volume;
        DataPro255{celli}{i,numCh} = canvasV;
        % Use the Mask on Fluorescence Images to process DIC Images
        for j = 1 : numCh-1
            cellthis = double(DataOri{i,celli}(:,:,j));
            sizeIm = size(cellthis);
            planeDICfit = fit([colN(:),rowN(:)],cellthis(:),'poly11','weights',box(:));
            planeDIC = planeDICfit.p00 + planeDICfit.p10*colN+planeDICfit.p01*rowN;
            normDIC = cellthis - planeDIC;
            Bgbox = box.*normDIC;
            Bgbox = Bgbox(:);
            Bgbox(Bgbox==0) = [];
            normDIC = (normDIC - mean(Bgbox)) / std(Bgbox);
            % padding
            muDIC = 0;
            stdDIC = 1;
            canvasDIC = normrnd(muDIC,stdDIC,[imSizeMax,imSizeMax]);
            [s1,s2] = size(cellthis);
            startRow = max(ceil((imSizeMax-1-s1)/2),1);
            startCol = max(ceil((imSizeMax-1-s2)/2),1);
            canvasDIC(startRow:startRow+s1-1,startCol:startCol+s2-1) = normDIC;
            canvasDIC255 = uint8(mat2gray(canvasDIC)*255);
            DataPro255{celli}{i,j} = canvasDIC255;
        end
    end
end

VolumeTime = nan(length(DataPro255),size(DataPro255{1},1));
for celli = 1:size(DataPro255,1)
    
    temp = cell2mat(DataPro255{celli,1}(:,2));
    
%
    VolumeTime(celli,1:length(temp)) = temp;% If two channel 4 to 3

end

save('VolumeTime','VolumeTime')
save('DataPro255','DataPro255','-v7.3')

