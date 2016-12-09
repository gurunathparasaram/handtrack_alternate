%ALTERNATE file for basic steps in HAND_TRACKING
%Author : PRG SWAMY
clc;
clear all;
close all;
clf; 

%read and display the image
img=imread('hand.jpg');
figure(1)
imshow(img);title('Original image');

%convert and display grayscale image
irgb=rgb2gray(img);
figure(2)
imshow(irgb);title('Grayscale image');

%use a 2D median filter to filter ou the salt and pepper noise
ifil=medfilt2(irgb);
figure(3)
imshow(ifil);title('2D median filtered image');

%perform thresholding on the image 
level=graythresh(ifil);
ibw=im2bw(ifil,level);
figure(4)
imshow(ibw);title('Thresholded image');

%Dilate and erode the image with diamond and square structural element of size 2,3 resply
sedil=strel('diamond',2);
idil=imdilate(ibw,sedil);
figure(5)
imshow(idil);title('Dilated image');
seerd=strel('square',3);
ierd=imerode(idil,seerd);
figure(6)
imshow(ierd);title('Eroded image');
img=ierd;

%Find the centroid of the hand using regionprops
statcent = regionprops(img,'centroid');
figure(7)
imshow(img); hold on;
for x = 1: numel(statcent)
    plot(statcent(x).Centroid(1),statcent(x).Centroid(2),'r.');
end
disp 'Centroid coordinates = '
disp statcent(1).Centroid(1);
disp statcent(1).Centroid(2);

%Find the orientation of the hand using regionprops
statorie=regionprops(img,'orientation');
figure(8)
imshow(img); hold on;
for x=1:numel(statorie)
    plot(statorie(x).Orientation);
end
disp 'Stat orientation :'
disp statorie.Orientation%Angle of orientation with respect to horizontal
%centx=statcent(1).Centroid(1);
%centy=statcent(1).Centroid(2);

%Find the angle of rotation(theta) 
if (0< statorie.Orientation(1) <90) %case of hand pointing outwards in 1st quadrant
      theta=(90-statorie.Orientation);%Angle to be rotated ACW for vertical alignment
elseif (-90 < statorie.Orientation(1) < 0)%case of hand pointing outwards in 4th quadrant
      theta=(-90-statorie.Orientation);%Angle to be rotated ACW for vertical alignment
end

%Rotate the hand ACW by theta to get alignment
imrot=imrotate(img,theta);
figure(9)
imshow(imrot);hold on; title('Rotated image');
statcentr = regionprops(imrot,'centroid');
for x = 1: numel(statcentr)
    plot(statcentr(x).Centroid(1),statcentr(x).Centroid(2),'r.');
end

centx=statcentr(x).Centroid(1);
centy=statcentr(x).Centroid(2);
BW = imrot;
%contour
figure(10)
BW2 = bwmorph(BW,'remove');
figure(11)
imshow(BW2)
% cont=imcontour(imrot,1);
[imht imwt]=size(BW2);
% % a=zeros(imht,imwt);
% img(round(cont(1,:)),round(cont(2,:))) = 1;
% for i=1:1:size(cont,2)
%     a(round(cont(2,i)),round(cont(1,i))) = 1;
% end
% figure(11)
% imshow(img);    
%!!!! Problem : contour image height is too big 

%to extract reference point by given algorithm
j=0;
tmp=[];
I=imrot;
for i=1:imht
    if(I(i,round(centy))==1)
        j=j+1;
        tmp(j)=i;
    end
    xr=max(tmp);
    yr=centy;
end 
figure(12)
imshow(imrot);title('Rotated image with reference point');hold on;
plot(yr,xr,'r.');
figure(13)
imshow(BW2);title('contour image with reference point ');hold on;
plot(yr,xr,'r.');

% centre = props.Centroid;
% [rows, cols] = find(BW2);
% dists = sqrt((rows-xr).^2 + (cols-yr).^2);
binaryImage = imrot;
fontSize = 13;
imshow(binaryImage);
axis on;
title('Binary Image', 'FontSize', fontSize);
title('With Boundaries, from bwboundaries()', 'FontSize', fontSize);
hold on;
boundaries = bwboundaries(binaryImage);
numberOfBoundaries = size(boundaries, 1); % Gives Number of Boundaries

for k = 1 : numberOfBoundaries
	thisBoundary = boundaries{k};
	plot(thisBoundary(:,2), thisBoundary(:,1), 'r', 'LineWidth', 3);
end
hold off;

% Define object boundaries
numberOfBoundaries = size(boundaries, 1)
message = sprintf('Found %d boundaries', numberOfBoundaries);
%uiwait(helpdlg(message));
boundary1 = boundaries{1};
boundary2 = [round(yr),round(xr)];
boundary1x = boundary1(:, 2);
boundary1y = boundary1(:, 1);
x1=1;
y1=1;
x2=1;
y2=1;
overallMaxDistance = inf; % Initialize.
% For every point in boundary 2, find the distance to every point in boundary 1.

for k = 1 : size(boundary1, 1)
	% Pick the next point on boundary 2.
	boundary2x = boundary2(1);
	boundary2y = boundary2(2);
	% For this point, compute distances from it to all points in boundary 1.
	allDistances = sqrt((boundary2x - boundary1x).^2 + (boundary2y - boundary1y).^2);
	% Find closest point, min distance.
	[maxDistance, indexOfMax1] = findpeaks(allDistances);
    windowSize = 6;
b = (1/windowSize)*ones(1,windowSize);
% b = [1 0.2];
a = [2];
maxDistance1 = medfilt1(allDistances,10);
% maxDistance1 = filter(b,a,allDistances);
% maxDistance1 = maxDistance1>20;
[maxDistance, indexOfMax] = findpeaks(maxDistance1);
% 	if maxDistance(k) < overallMaxDistance
		x1 = boundary1x(indexOfMax);
		y1 = boundary1y(indexOfMax);
		x2 = boundary2x;
		y2 = boundary2y;
		overallMaxDistance = maxDistance;
	end
% end
% Find the overall min distance
% maxDistance = max(maxDistance);
% Report to command window.
% maxDistance

% Draw a line between point 1 and 2
for i = 1:size(x1,1)
line([x1(i), x2], [y1(i), y2], 'Color', 'b', 'LineWidth', 3);
end

% [pks,locs] = findpeaks(allDistances);
figure,
findpeaks(maxDistance)
text(indexOfMax+.02,maxDistance,num2str((1:numel(maxDistance))'))

