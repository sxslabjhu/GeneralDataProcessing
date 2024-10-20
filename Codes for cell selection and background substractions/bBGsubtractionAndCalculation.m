

 clear
 clc
 close all

dataOutput = load('rawSingleCells.mat');
dataOutput = dataOutput.dataOutput;

%Define Fucci Channels.
fucciChannels = [1];

fucciOutput = {nan(size(dataOutput)),nan(size(dataOutput))};% cell=A{[],[]} put two matrix into one cell
Intensity = fucciOutput;
fucciOutput2 = fucciOutput;

for j = 1:length(dataOutput(:,1))%j is the time 
    for k = 1:length(dataOutput(1,:))% k is the cell number
        if ~isempty(dataOutput{j,k})
            for i = 1:length(fucciChannels)% i the single number
                %complete image processing routine for each fucci image
                [fucciOutput{i}(j,k), fucciOutput2{i}(j,k)]= intensityFunction(dataOutput{j,k}(:,:,fucciChannels(i)));
            end
        end
    end
end
%save('fucciOutput','fucciOutput')
 meanIntensity=fucciOutput{1,1};%% For arcklight
 meanIntensity=meanIntensity';
 %Intensity(:,6)=[];
 
 Intensity=fucciOutput2{1,1};%% For arcklight
 Intensity=Intensity';
 save('Intensity','Intensity','meanIntensity')

% 
% figure
% plot(Intensity'./Intensity(:,1)')

figure
%temp = Intensity'./Intensity(:,1)';
temp = meanIntensity'./meanIntensity(:,1)';
plot(temp)
hold on
%errorbar(mean(temp,2,'omitnan'),std(temp,[],2,'omitnan'),'linewidth',2,'color','k')

function [meanIntensity,sumI] = intensityFunction(rawROI)

%Fit a plane (linear function) to background. Subtract plane, and sum.
d = 10;%size of boundary box
[colN,rowN] = meshgrid(linspace(1,length(rawROI(1,:)),length(rawROI(1,:))),...
    linspace(1,length(rawROI(:,1)),length(rawROI(:,1))));
box = logical(gt(colN,max(colN(:))-d)+lt(colN,d+1)+gt(rowN,max(rowN(:))-d)+lt(rowN,d+1));%gt(A,B) is an alternate way to execute A > B.
%colN(:) make whole matrix a coloum vector,it(A,B) is an alternate way to execute A < B.
fitSurface = fit([colN(:),rowN(:)],rawROI(:),'poly11','weights',box(:));
fitPlane = fitSurface.p00 + fitSurface.p10*colN+fitSurface.p01*rowN;

normROI = rawROI-fitPlane;

imtemp = imgaussfilt(normROI);
BgAll = mean(imtemp(:));
StdAll = std(imtemp(:));

MaskAllDNA = imtemp>BgAll + StdAll;
MaskAllDNA = bwareaopen(MaskAllDNA,200);
MaskAllDNA = imfill(MaskAllDNA,'holes');
%MaskAllDNA = imclearborder(MaskAllDNA);


meanIntensity = mean(normROI(MaskAllDNA));
sumI = sum(normROI(:));

end




