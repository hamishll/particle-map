close all;
clc;
clear;
%%// Input image
IMG = imread('img_batch_02/comp.png'); %%// Available in the MATLAB image library
figure,imshow(IMG);
[imgy, imgx] = size(IMG);

%%// Conversion to BW image through thresholding
WB = im2bw(IMG, 0.67);
BW = imcomplement(WB);
figure,imshow(BW);

%%// Get all measurements into one structure
s = regionprops(BW, 'Area', 'Perimeter', 'EquivDiameter', 'Centroid', 'BoundingBox');
%%// We can then pull each component as ARRAYS (instead of STRUCTURES)
example = [s.Area];

%%// discard particles with area less than N pixels
s2 = struct('Area',{}, 'Perimeter',{}, 'EquivDiameter',{}, 'Centroid',{}, 'BoundingBox',{});
for i=1:length([s.Area])
    if s(i).Area > 5
       s2 = [s2, s(i)];
    else
    end
end

%%// THIS PART IS REDUNDANT
%%// Get the scalar distances around the boundaries of the regions/blobs
P  = regionprops(BW, 'Perimeter');
%%// Get scalar values that specifies the diameter of a circle with the same
%%// area as the regions/blobs
D  = regionprops(BW, 'EquivDiameter');
A  = regionprops(BW, 'Area');
%%// THIS PART IS REDUNDANT

%%// Getting a distribution of areas [*]
%allDiameters = [D.EquivDiameter];
allAreas = [s2.Area];
bins = 500;
%[diamDistribution, binDiameters] = hist(allDiameters, bins);
[areaDistribution, binAreas] = hist(allAreas, bins);
figure,bar(binAreas, areaDistribution, 'BarWidth', 1.0);
axis([0,inf,0,inf]);

%%// Finding the centre of mass for each particle
COM = regionprops(BW, 'Centroid');
%%// Also redundant now ^^^

x_centroid = [0,0,0,0];
y_centroid = [0,0,0,0];

for i=1:length([s2.Area])
    x_centroid(i) = s2(i).Centroid(1);
    y_centroid(i) = s2(i).Centroid(2);
    size_centroid(i) = s2(i).Area;
end
% sqrt(A/pi)
% we need to invert y_centroid to display a scatter that isn't a reflection
for i=1:length(y_centroid);
    y_centroid_inv(i) = imgy - y_centroid(i);
end   

figure,scatter(x_centroid, y_centroid_inv, (size_centroid/2));
axis([0,imgx,0,imgy]);

%%// X-Distribution
binsPos = 20;
[posDistribution, binDiameters] = hist(x_centroid, binsPos); 
figure,bar(binDiameters, posDistribution, 'BarWidth', 1.0);
axis([0,imgx,0,inf]);

%%// Y-Distribution



%%// Getting a distribution of particle size across x

%%// Let's see if we can get a log-normal fitting

%%// Note that you can do this:    
%%// s = regionprops(L, 'Area', 'Centroid', 'BoundingBox'); 




%%// Get list of pixels based on their labeling. 
%%// Basically the indices of the structs produced by regionprops refer to the labels.
pixel_list = regionprops(BW, 'PixelIdxList');

%%// Let us find out information about the Nth blob (blob that is labeled as N)
N = 1;

%%// 1. List of pixel coordinates as linear indices
blob1_pixel_list = pixel_list(N).PixelIdxList;
% this line ^ tells us a lot about syntax 

%%// Create an image of the same size as the original one and showing the
%%blob labeled as 1
blob1 = false(size(BW));
blob1(blob1_pixel_list) = true;
figure,imshow(blob1);

%%// Perimter of blob -1 
blob1_perimeter = P(N).Perimeter;

%%// Equivalent diameter of blob -1 
blob1_equivdiameter_values = D(N).EquivDiameter;

%%// Get perimeter and diameter values for all the blobs
perimeter_values = struct2array(P);
diameter_values = struct2array(D);