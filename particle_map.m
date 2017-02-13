close all;
clc;
clear;

%%// TO-DO
% noise reduction within MATLAB
% multiple threshold sampling
% view intensity histogram, prior to BW

%%// FIGURE(1): Input image
IMG = imread('img_batch_02/comp_nr.png'); %%// Available in the MATLAB image library
figure(1); subplot(1,3,1);
imshow(IMG);
[imgy, imgx] = size(IMG);
b = [1,2,3];
%%/////////////////////////////////////////////////////////
%%//  Viewing histogram of intensities
%%/////////////////////////////////////////////////////////
% need to create a 1D array from IMG
allIntensities = double(IMG(:)');   %IMG(:)';
binsInt = 254;
[DistributionInt, binInt] = hist(allIntensities, binsInt);
figure(2); subplot(4,1,1);
bar(binInt, DistributionInt, 'BarWidth', 1.0);
axis([0,254,0,inf]);
findpeaks(DistributionInt, 1, 'MinPeakDistance', 20);

%%/////////////////////////////////////////////////////////
%%// Converting image to binary; collecting data
%%/////////////////////////////////////////////////////////

%%// FIGURE(2): Conversion to BW image through thresholding
%WB = im2bw(IMG, 0.67);
%BW = imcomplement(WB);
%figure,imshow(BW);

%%// Alternative FIGURE(2): thresholding between two values
BW = IMG;
th1 = 0; th2 = 166;
th_range = (BW > th1 & BW <= th2);
BW(th_range) = 1;
BW(~th_range) = 0;
BW = logical(BW);
%figure,imshow(BW)

%%// Get all measurements into one structure
s = regionprops(BW, 'Area', 'Perimeter', 'EquivDiameter', 'Centroid', 'BoundingBox');
%%// We can then pull each component as ARRAYS (instead of STRUCTURES)
example = [s.Area];

%%/////////////////////////////////////////////////////////
%%// Filtering out particles with area less than N pixels
%%/////////////////////////////////////////////////////////

s2 = struct('Area',{}, 'Perimeter',{}, 'EquivDiameter',{}, 'Centroid',{}, 'BoundingBox',{});
for i=1:length([s.Area])
    if s(i).Area > 10
       s2 = [s2, s(i)];
    else
    end
end

figure(1); subplot(1,3,2);
imshow(BW)
hold on
centroids = cat(1, s2.Centroid);
plot(centroids(:,1),centroids(:,2), 'r.')
hold off;

%%/////////////////////////////////////////////////////////
%%//  Size distribution
%%/////////////////////////////////////////////////////////

%FIGURE(3): Getting a distribution of sizes:
allAreas = [s2.Area];
bins = 2000;
[areaDistribution, binAreas] = hist(allAreas, bins);
figure(2); subplot(4,1,2);
bar(binAreas, areaDistribution, 'BarWidth', 1.0);
axis([0,1000,0,inf]);

%%/////////////////////////////////////////////////////////
%%// Plotting scatters
%%/////////////////////////////////////////////////////////

for i=1:length([s2.Area])
    x_centroid(i) = s2(i).Centroid(1);
    y_centroid(i) = s2(i).Centroid(2);
    size_centroid(i) = s2(i).Area;
end

% FIGURE(4): x,y scatter
% we need to invert y_centroid to display a scatter that isn't a reflection
for i=1:length(y_centroid);
    y_centroid_inv(i) = imgy - y_centroid(i);
end   
figure(2); subplot(4,1,3);
scatter(x_centroid, y_centroid_inv, (size_centroid/8));
axis([0,imgx,0,imgy]);

%%/////////////////////////////////////////////////////////
%%// Plotting distributions
%%/////////////////////////////////////////////////////////

% FIGURE(5): X-Distribution
binsPos = 20;
[posDistribution, binDiameters] = hist(x_centroid, binsPos); 
figure(2); subplot(4,1,4);
bar(binDiameters, posDistribution, 'BarWidth', 1.0);
axis([0,imgx,0,inf]);

%%// Y-Distribution

%%// Getting a distribution of particle size across x

%%// Let's see if we can get a log-normal fitting



%%/////////////////////////////////////////////////////////
%%// STORING AND RECALLING SPECIFIC BLOBS
%%/////////////////////////////////////////////////////////

%%// Get list of pixels based on their labeling. 
%%// The indices of the structs produced by regionprops refer to the labels
pixel_list = regionprops(BW, 'PixelIdxList');

%%// Let us find out information about the Nth blob (blob that is labeled as N)
N = 80;

%%// 1. List of pixel coordinates as linear indices
blobN_pixel_list = pixel_list(N).PixelIdxList;
% this line ^ tells us a lot about syntax 

%%// FIGURE(6): Create an image of the same size as the original one and showing the
% blob labeled N
blobN = false(size(BW));
blobN(blobN_pixel_list) = true;
figure(1); subplot(3,1,3);
imshow(blobN);

% Perimter of blob N
blobN_perimeter = s2(N).Perimeter;
% Equivalent diameter of blob N
blobN_equivdiameter_values = s2(N).EquivDiameter;