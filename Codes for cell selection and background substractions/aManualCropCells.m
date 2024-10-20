 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is associated with the paper titled "Cytoskeletal activation of
% NHE1 regulates mechanosensitive cell volume adaptation and proliferation"
% by Ni, et al. 
% The code is can be used for academic and research purposes only.
% For inquiries, please contact Qin Ni or Sean X. Sun at Johns Hopkins
% University. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cutting Cells From tiff files, manual

clear
clc
close all

%Define Number of channels. This script assumes the DIC channel is always
%first in the output order.
numChannel = 1;

%Identify images
rector = dir('*.tif');
numTimepoints = length(rector)/numChannel;

%Initialize output file
dataOutput = {};%each row corresponds to a timepoint, each column to a cell

%Save first image to select cells
DIC1 = imread(rector(1).name);

DIC2 = imread(rector(2).name);
%Initialize while loop to select cells
stoppor = 0;
cellIndex = 0;
a = length(rector(:,1))/2;
p = 0;


while ~stoppor
    cellIndex = cellIndex + 1;
    imagesc(DIC1)
    colorbar
    %colormap gray

    [x,y] = ginput(2); %select cell by clicking corners of rectangle,
                       %finish loop (no more good cells) by double-clicking
    
    x1 = floor(min(x));
    x2 = ceil(max(x));
    y1 = floor(min(y));
    y2 = ceil(max(y));
    
    if x2-x1 < 5 %meaning you have double-clicked, no more good cells
        stoppor = 1; %script ends
        
        
    elseif p > a
        stoppor = 1;
        
    else
        DIC1(y1:y2,x1:x2) = mean(DIC1(:));%remove selected cell so you don't select it twice
        currCell = {}; %initialize selected cell
        
        %check that the boundary is good for every timepoint and save the
        %data
        i = 0;  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        stoppor2 = 0;
        while lt(i,numTimepoints)&&(~stoppor2)% %%for scalar calculation, when a%b, only when a=1 b=1 result=1; % 1t(A B) is equal to A<B.
            i = i + 1;
            p = p + 1;
            timepointBlock = cell(numChannel,1);%saves all the channels from the current timepoint
            for j = 1:numChannel
               timepointBlock{j} = double(imread(rector((i-1)*numChannel+j).name));
            end
            
            %check cell is in ROI
            imagesc(timepointBlock{1}(y1:y2,x1:x2))
           title(num2str(i))
            %colormap gray
            
            %press 'enter' or any other key to keep the box and '0' to reset the box
            waitforbuttonpress;
            keepBox = get(gcf,'CurrentCharacter');% Record the key you press in the keyboard
            if keepBox=='0'                
                checkROI = 0;
                while ~checkROI
                    %show original box and choose new box
                    imagesc(timepointBlock{1})
                    colorbar
                    %colormap gray
                    rectangle('Position',[x1,y1,x2-x1,y2-y1])
                    title(num2str(i))

                    [x,y] = ginput(2); %select cell by clicking corners of rectangle,
                    %finish loop (no more good cells) by double-clicking
                    
                    x1 = floor(min(x));
                    x2 = ceil(max(x));
                    y1 = floor(min(y));
                    y2 = ceil(max(y));
        

                    if x2-x1 < 5 %meaning you have double-clicked,
                        checkROI = 1;
                        stoppor2 = 1; %stop following this cell
                        
                    else
                        
                        imagesc(timepointBlock{1}(y1:y2,x1:x2))
                        title(num2str(i))
                        %colormap gray
                        waitforbuttonpress;
                        checkROI = get(gcf,'CurrentCharacter');
                        checkROI = ~isequal(checkROI,'0');% is checkROI=0, isequal is 1,~isequal is 0;
                    end
                end
            end
            if ~stoppor2
                currBlock = nan(y2-y1+1,x2-x1+1,numChannel); 
                for k = 1:numChannel
                    currBlock(:,:,k) = timepointBlock{k}(y1:y2,x1:x2);  
                end 
                currCell = cat(1,currCell,currBlock); 
            end   
        end 
        dataOutput(1:length(currCell),cellIndex) = currCell; %#ok<SAGROW>    
    end          
end         
close all      


%% save this cell
savefilename = dir('../Combine/rawSingleCells*');

if isempty(dataOutput) 
    disp('Error: no file!')
else
    if isempty(savefilename) 
        nextname = 1
    else
        allname = zeros(length(savefilename),1);
        for ni = 1:length(savefilename)
            nametemp = savefilename(ni).name;
            endnamei = find(nametemp == '.',1) - 1;
            allname(ni) = str2num(nametemp(15:endnamei));
        end
        
        nextname = max(allname) + 1
    end
    
    save(['../Combine/rawSingleCells',num2str(nextname)],'dataOutput')
    
    clear dataOutput
    
end



