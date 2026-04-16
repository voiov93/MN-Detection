%% Vivian's MN finding script


%% STEP 1: Load image

clc
clear all
close all

%I_orig = imread("Z:\2025_ALL_USER_IMAGES_HERE\Vivian Cho\i880\gsk_12_hr\crop_test.tif");
%I_orig = imread("C:\Users\vivia\Downloads\MAX_40x_zstack_tilescan.tif");
%I_orig = imread("Z:\2025_ALL_USER_IMAGES_HERE\Vivian Cho\i880\gsk_12_hr\AVG_40x_zstack_tilescan.tif");
%I_orig = imread("Z:\2025_ALL_USER_IMAGES_HERE\Vivian Cho\i880\gsk_12_hr\Combined.tif");

I_orig = imread("Z:\2026_ALL_USER_IMAGES_HERE\Vivian Cho\4_15\dish1\MAX_C2-dish1_40x_1_stitch.tif");

% time series image MCF10As with GFP
%I_orig = imread("Z:\2025_ALL_USER_IMAGES_HERE\Vivian Cho\i880\gsk_16_hr_10A\Time Slices\time_tile_40x_10015.tif");
I = im2gray(I_orig); 
I = imlocalbrighten(I); 
subplot(3,3,1);
imshow(I)
title('imlocalbrighten')
hold on

%% STEP 2: Subtract background noise based on histogram and Smooth

% histogram: 
h = histogram(I, 10); 
counts = h.Values; 
edges = h.BinEdges; 
x = find(I<edges(3));
I(x) = 0;

img_blur = imgaussfilt(I, 2); 
imshow(img_blur, [])
I = img_blur; 
subplot(3,3,2); 
imshow(I);
title('subtract background')
hold on
%% STEP 3: MULTITHRESH
thresh = multithresh(I,20);

labels = imquantize(I,thresh); 
I_labelsRGB = label2rgb(labels); 
%thresh_labels = labels>13; 
thresh_labels = labels>0; 
j = labels .* thresh_labels;

% Get rid of the bottom 13 thresholds. 
I_20 = label2rgb(j); 
%imshow(I_20);
subplot(3,3,3); 
imshow(I_labelsRGB); 
title('multithresh 20')
hold on

subplot(3,3,4);
imshow(I_20);
title('Only top 7 threshold')
hold on

%% STEP 4:  Binarize
subplot(3,3,5)
I_g = im2gray(I_20);
BW = imbinarize(I_g);
subplot(3,3,5)
imshow(BW);
title('binarize')
hold on


%% edge detector pt 2: 
% % SKIP THIS 
% [Gx,Gy] = imgradientxy(I);
% [Gmag,Gdir] = imgradient(Gx,Gy);
% 
% Gmag_t = Gmag>300; 
% Gmag_two = Gmag.*Gmag_t; 
% imshowpair(Gmag,Gmag_two,'montage');

%% STEP 5  % Edge detector and defining Previous Threshold features.  
I = BW; 

[~, threshold] = edge(I, 'Sobel'); 
fudgeFactor = 0.5; 
% sobel is an EDGE detector based on gradient
BWs = edge(I, 'sobel', threshold*fudgeFactor); 

% dilate the sobel image. 
se90 = strel('line', 3,90); 
se0 = strel('line', 3,0); 

BWsdil = imdilate(BWs, [se90, se0]); 
%figure(4); 
%imshow(BWsdil); 

% Fill interior Gaps
BWdfill = imfill(BWsdil, 'holes');

% erode the sobel image
%se_er = strel("line", 3,70);
BWsero = imerode(BWdfill, [se0, se90]);
% clean up edges
J = imclearborder(BWsero); 
subplot(3,3,6)
imshow(J); 
title('edge detector')
hold on

%% STEP 6: ONLY keep nuclei and eliminate mitotic catastrophes

CC = bwconncomp(J); 

[rows, columns] = size(J);  % size of CC is 1x1
 

% to filter for large nuclei. 
stats = regionprops("table",CC,"Centroid", ...
    "MajorAxisLength","MinorAxisLength", "Perimeter"); 
diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);

num_bins = 20;
h = histogram(diameters,num_bins);
counts = h.Values; 
edges = h.BinEdges; 
%x = find(diameters<edges(3));
%diameters(x) = 0;
%close all

% uses otsuthresh to find the bin number separating nuclei from MN.
%T = floor(otsuthresh(diameters)*num_bins)-2; 
T = 6; % edge
% stores nucli that are within the histogram bin filter
% NOT to filter for micronuclei. Just a very general way of filtering
nuclei_filt1 = find(diameters > edges(T) & diameters < max(diameters)+1); 

% this code filters out two nuclei that are too close
nuclei_not_too_close = [];
nuclei_too_close = [];
for i = 1:length(nuclei_filt1)
    n_i = nuclei_filt1(i); % holds the nuclei object number 
    center = stats.Centroid(n_i,:);
    % distance between this nuclei (n_i) and nuclei_filt_1
    dx = center(1,1) - stats.Centroid(nuclei_filt1,1); 
    dy = center(1,2) - stats.Centroid(nuclei_filt1,2); 
    distance = sqrt(dx.^2+dy.^2);
    
    % too_close is index position in distance. 
    too_close = find(distance<diameters(n_i)*0.9 & distance~=0);
    if isempty(too_close) 
        %reassign
        nuclei_not_too_close = [nuclei_not_too_close, n_i]; 
    else 
        too_close_nuclei = nuclei_filt1(too_close); % find nuclei object #
        both_too_close = [n_i too_close_nuclei']; % BOTH nuclei and nuclei close to it stored
        nuclei_too_close = [nuclei_too_close, both_too_close];
    end
end
% stores nuclei that are not too close to each other
nuclei_filt2 = setxor(nuclei_too_close, nuclei_filt1);
%imshow(cc2bw(CC, ObjectsToKeep = c))


% filter against less circular nulcei. 
stats_c = regionprops("table", CC, "Circularity"); 
circular_general = find(stats_c.Circularity > 0.6);

% contains only normal shaped nuclei, no abnormal ones. 
nuclei_idx = intersect(nuclei_filt2, circular_general);
subplot(3,3,7)
imshow(cc2bw(CC, ObjectsToKeep = nuclei_idx)); 
title('ONLY nuclei')
hold on


%% STEP 7: Only keeping Micronuclei by calculating distance
% This is different than using a histogram in the previous step. 

% separate POSSIBLE MICRONUCLEI (REFERENCED AS "DOTS") and nuclei 
centroids = stats.Centroid; 
nuclei_centroids = centroids(nuclei_idx,:);

% create an empty cell to store shapes that are not nuclei
dot_idx_cell = cell(length(centroids)-length(nuclei_idx),1);
for i = 1:length(centroids)
    if i ~= nuclei_idx
        dot_idx_cell{i,1} = i;
    end
end

% convert cell to array only keep nuclei below histogram edge threshold
dot_idx_pre_filt1 = cell2mat(dot_idx_cell); 
dot_idx2 = dot_idx_pre_filt1(find(diameters(dot_idx_pre_filt1) < edges(T))); % find the index of the vector (diam...) 
% FILTER BASED ON CIRCULARITY
dot_idx = dot_idx2(find(stats_c.Circularity(dot_idx2)>0.8)); 
dot_centroids = centroids(dot_idx,:);


% filter 1
% mn_idx should only contain points from dots_idx that ARE MN. 
mn_idx = cell(length(nuclei_idx),1);

% this will store nuclei that ARE NOT mitotic catastrophe
nuclei_idx_real = nuclei_idx;

for i = 1:length(nuclei_idx)
    % find this iteration's nuclei center 
    temp_nuclei_centroid = nuclei_centroids(i,:);
    % find this iteration's nuclei diameter
    temp_nuclei_diameter = diameters(nuclei_idx(i));
    
    % (1) dot centroid is within a certain radius of nuclei
    dx = dot_centroids(:,1) - temp_nuclei_centroid(1,1); 
    dy = dot_centroids(:,2) - temp_nuclei_centroid(1,2); 
    temp_distance = sqrt(dx.^2+dy.^2);
    
    % if shape is within a X radius from nuclei center, it is mn 
    %constraint = temp_nuclei_diameter/2*1.65;  % distance constraint
    constraint = temp_nuclei_diameter/2*2.5;
    dot_idx_idx =   abs(temp_distance) < constraint; 
    mn_indexes = dot_idx(dot_idx_idx);
    mn_idx{i, 1} = mn_indexes; 

    % (2) if a cell has more than 3 MN, it is mitotic catastrophe
    if length(mn_idx{i,:})>3
        mn_idx{i,:} = 0;
        nuclei_idx_real(i) = 0;
    end
end

% convert this cell into array
%mn_idx = mn_idx(~cellfun('isempty',mn_idx));
idx_not_zero = ~cellfun(@(x) any(x == 0), mn_idx);
mn_idx = mn_idx(idx_not_zero);
mn_idx_array = cell2mat(mn_idx);
mn_idx_array = mn_idx_array(mn_idx_array>0);


nuclei_idx_real = nuclei_idx_real(nuclei_idx_real>0);
% display only the MN 
subplot(3,3,8)
imshow(cc2bw(CC, ObjectsToKeep = mn_idx_array));
title('ONLY MN')
hold off

%% MN and Nuclei last Step and %% plot over original image

three = [mn_idx_array; nuclei_idx_real];

%imshow(cc2bw(CC, ObjectsToKeep=three));


figure(2); 
imshow(I_orig); 
hold on; 
plot(centroids(mn_idx_array,1), centroids(mn_idx_array,2), '*m');
hold on; 
%viscircles(centroids(mn_idx_array,:), diameters(mn_idx_array)/2, 'LineWidth', 0.5, 'Color', 'm');
%hold on;
%plot(centroids(nuclei_idx_real,1), centroids(nuclei_idx_real',2), '*g');
%hold on; 
plot(centroids(nuclei_idx_real,1), centroids(nuclei_idx_real',2), '*g');
hold on; 
viscircles(centroids(nuclei_idx_real,:), diameters(nuclei_idx_real)/2, 'LineWidth', 0.5, 'Color', 'g');
%hold on; 
title('MN and Nuclei on Original Image - no corrections')
hold off;
%size(mn_idx_array)
%size(nuclei_idx_real)


%% MN and Nuclei Corrections based on original image. 
% some mn don't have the right intensity compared to the nucleus

% find the centroids, diameters, pixels of the nuclei from last step
cc_centroids = centroids(nuclei_idx_real);
cc_diameters = diameters(nuclei_idx_real);
cc_nu_pxid = CC.PixelIdxList(nuclei_idx_real);

% mn_idx is the cell housing the mn from the last step: REFERENCE IT

% this is a empty cell to house all mn I want to KEEP 
keeper_mn_idx = cell(length(mn_idx_array),1); 


% loops in respect to the nuclei 
for i = 1:length(nuclei_idx_real)
    % this current iteration's nuclei pixels
    nu_pxid_list = cc_nu_pxid{:,i}; 
    % find average nuclei intensity based on original pixels
    average_intensity_nu = mean(I_orig(nu_pxid_list)); 

    % holds the micronuclei belonging to current interation's nuclei 
    mn_temp = mn_idx{i,:};
    y = isempty(mn_temp); % if this is 1, that means it is empty

    if isempty(mn_temp) < 1 & ismember(mn_temp, 0) < 1
        cc_mn_pxid = CC.PixelIdxList(mn_temp);
        %mn_vals = [];
        x = mn_temp;
        for j = 1:length(mn_temp)
            mn_temp_pxid = cc_mn_pxid{:,j}; 
            int_vals_mn = mean(I_orig (mn_temp_pxid));
            % if the MN has a lower intensity than nuclei, it is deleted
            % 0.89 works on the smaller one 
            if int_vals_mn < average_intensity_nu*0.63
                % if less than 50% of nucleus intensity, eliminate
                x(ismember(x, mn_temp(j))) = [];
                keeper_mn_idx{i, :} = x; 
            else 
                % if not, keep the micronuclei. 
                keeper_mn_idx{i, :} = x;
            end
        end  
    end

end

% pad short rows: 
maxlength = max(cellfun(@length, keeper_mn_idx));
keeper_mn_idx_padded = cellfun(@(x) [x; NaN(maxlength - length(x),1)], keeper_mn_idx, 'UniformOutput', false);

keeper_mn_idx_array = cell2mat(keeper_mn_idx_padded);
keeper_mn_idx_array = keeper_mn_idx_array(keeper_mn_idx_array>0);

% plot new MN over the image: 
figure(3); 
imshow(I_orig); 
hold on; 
plot(centroids(keeper_mn_idx_array,1), centroids(keeper_mn_idx_array,2), '*m');
hold on; 
%plot(centroids(nuclei_idx_real,1), centroids(nuclei_idx_real',2), '*g');
%hold on; 
plot(centroids(nuclei_idx_real,1), centroids(nuclei_idx_real',2), '*g');
hold on; 
viscircles(centroids(nuclei_idx_real,:), diameters(nuclei_idx_real)/2, 'LineWidth', 0.5, 'Color', 'g');
hold on; 
title('MN and Nuclei on Original Image - corrected');



% mn_idx_array may have repeats: 
mn_idx_array = unique(mn_idx_array);
size(mn_idx_array)
size(nuclei_idx_real)

% compare to the intensity pictures
figure(8)
subplot(1,2,1)
imshow(I_labelsRGB)
subplot(1,2,2)
imshow(I_20)