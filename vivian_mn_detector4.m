%% Vivian MN detector - cleaned integrated version
% RULES:
% 1) Nuclei are detected from H1
% 2) Standard MN candidates are detected from H1
% 3) Extra MN candidates can also come from bright, round Emerin circles
% 4) Emerin-circle candidates must overlap H1 somewhere
% 5) Each MN is assigned to the nearest nucleus
% 6) If MN centroid is inside assigned parent nucleus, reject it
% 7) If a nucleus has 3 or more MN, remove that nucleus and all its MN
% 8) Exclude nuclei and MN near the image border
%
% Display:
% - Overlay shown on Emerin image
% - Green = valid nuclei
% - Red = final MN
% - Yellow = candidate MN before final exclusion

clc
clear
close all

%% ---------------- LOAD IMAGES ----------------
% PICK ONE SET OF PATHS ONLY

H1 = imread('Z:\2026_ALL_USER_IMAGES_HERE\Vivian Cho\4_15\dish4\1\MAX_C2-dish4_40x_1_stitch.tif');
Emerin = imread('Z:\2026_ALL_USER_IMAGES_HERE\Vivian Cho\4_15\dish4\1\MAX_C1-dish4_40x_1_stitch.tif');

% If you want the "_3" images instead, comment the two lines above and use:
% H1 = imread('Z:\2026_ALL_USER_IMAGES_HERE\Vivian Cho\4_15\dish1\3\MAX_C2-dish1_40x_3_stitch.tif');
% Emerin = imread('Z:\2026_ALL_USER_IMAGES_HERE\Vivian Cho\4_15\dish1\3\MAX_C1-dish1_40x_3_stitch.tif');

if ndims(H1) == 3
    H1 = rgb2gray(H1);
end
if ndims(Emerin) == 3
    Emerin = rgb2gray(Emerin);
end

H1 = mat2gray(H1);
Emerin = mat2gray(Emerin);

%% ---------------- PARAMETERS ----------------
% Standard MN thresholds from H1
minMicronucleusArea = 5;
maxMicronucleusArea = 150;
minSolidityMN = 0.60;
minCircularityMN = 0.25;

% Nucleus thresholds from H1
minNucleusArea = 200;
maxNucleusArea = inf;
minSolidityNucleus = 0.75;
minCircularityNucleus = 0.20;

% Extra bright round Emerin-circle candidate settings
useBrightRoundEmerinCandidates = true;
minBrightEmerinArea = 20;
maxBrightEmerinArea = 250;
minBrightEmerinCircularity = 0.60;
minBrightEmerinSolidity = 0.80;
minBrightEmerinMeanIntensity = 0.20;

% Rescue very bright intranuclear emerin-positive MN-like objects
allowVeryBrightIntranuclearEmerinMN = true;
minRescueEmerinMeanIntensity = 0.35;   % tune
minRescueEmerinCircularity = 0.70;     % tune
minRescueEmerinSolidity = 0.80;        % tune
minRescueEmerinArea = 10;              % tune
maxRescueEmerinArea = 200;             % tune

% Preprocessing
gaussianSigma = 1.0;
adaptiveSensitivityH1 = 0.4;
adaptiveSensitivityEmerin = 0.3;

% Emerin background subtraction
bgRadius = 20;

% Parent assignment
maxDistanceToNearestNucleus = 60;   % pixels

% Border exclusion
borderMargin = 5;

% Reject MN if inside assigned parent nucleus
rejectMNInsideAnyNucleus = true;

% Deduplicate candidates closer than this distance
candidateMergeDistance = 5;   % pixels

% Watershed
doWatershedH1 = true;
doWatershedH1Candidates = false;

% Display
showDebugFigures = true;

%% ---------------- PREPROCESS ----------------
H1_smooth = imgaussfilt(H1, gaussianSigma);
H1_enh = adapthisteq(H1_smooth);

Emerin_smooth = imgaussfilt(Emerin, gaussianSigma);
se_bg = strel('disk', bgRadius);
Emerin_sub = imtophat(Emerin_smooth, se_bg);
Emerin_enh = adapthisteq(mat2gray(Emerin_sub));

if showDebugFigures
    figure; imshow(H1, []); title('Original H1');
    figure; imshow(H1_enh, []); title('Enhanced H1');
    figure; imshow(Emerin, []); title('Original Emerin');
    figure; imshow(Emerin_sub, []); title('Emerin background-subtracted');
    figure; imshow(Emerin_enh, []); title('Enhanced Emerin');
end

%% ---------------- SEGMENT H1 FOR NUCLEI ----------------
BW_H1_nuc = imbinarize(H1_enh, 'adaptive', ...
    'ForegroundPolarity', 'bright', ...
    'Sensitivity', adaptiveSensitivityH1);

BW_H1_nuc = bwareaopen(BW_H1_nuc, minMicronucleusArea);
BW_H1_nuc = imfill(BW_H1_nuc, 'holes');

if doWatershedH1
    D = -bwdist(~BW_H1_nuc);
    D(~BW_H1_nuc) = -Inf;
    D2 = imhmin(D, 1);
    L = watershed(D2);
    BW_H1_nuc(L == 0) = 0;
end

if showDebugFigures
    figure; imshow(BW_H1_nuc); title('H1 nuclei segmentation');
end

%% ---------------- SEGMENT H1 FOR MN CANDIDATES ----------------
BW_H1_cand = imbinarize(H1_enh, 'adaptive', ...
    'ForegroundPolarity', 'bright', ...
    'Sensitivity', adaptiveSensitivityH1);

BW_H1_cand = bwareaopen(BW_H1_cand, minMicronucleusArea);
BW_H1_cand = imfill(BW_H1_cand, 'holes');

if doWatershedH1Candidates
    D = -bwdist(~BW_H1_cand);
    D(~BW_H1_cand) = -Inf;
    D2 = imhmin(D, 1);
    L = watershed(D2);
    BW_H1_cand(L == 0) = 0;
end

%% ---------------- SEGMENT EMERIN SUPPORT MASK ----------------
BW_E_support = imbinarize(Emerin_enh, 'adaptive', ...
    'ForegroundPolarity', 'bright', ...
    'Sensitivity', adaptiveSensitivityEmerin);

BW_E_support = bwareaopen(BW_E_support, minMicronucleusArea);
BW_E_support = imfill(BW_E_support, 'holes');

% Break weak bridges a little
se_break = strel('disk', 1);
BW_E_support = imerode(BW_E_support, se_break);
BW_E_support = imdilate(BW_E_support, se_break);

% Watershed to split touching objects
D = -bwdist(~BW_E_support);
D(~BW_E_support) = -Inf;
D2 = imhmin(D, 1);   % try 1, 2, or 3
L = watershed(D2);
BW_E_support(L == 0) = 0;

if showDebugFigures
    figure; imshow(BW_E_support); title('Emerin support mask after watershed');
end

%% ---------------- REGIONPROPS ----------------
statsH1 = regionprops(BW_H1_nuc, H1_enh, ...
    'Area', 'Perimeter', 'Solidity', 'Centroid', 'PixelIdxList', 'BoundingBox');

statsMN = regionprops(BW_H1_cand, H1_enh, ...
    'Area', 'Perimeter', 'Solidity', 'Centroid', 'PixelIdxList', 'BoundingBox');

statsE = regionprops(BW_E_support, Emerin_enh, ...
    'Area', 'Perimeter', 'Solidity', 'Centroid', 'PixelIdxList', 'BoundingBox', 'MeanIntensity');

imgH = size(H1,1);
imgW = size(H1,2);

%% ---------------- CLASSIFY NUCLEI FROM H1 ----------------
nH1 = numel(statsH1);

areasH1 = [statsH1.Area]';
perimH1 = [statsH1.Perimeter]';
solH1 = [statsH1.Solidity]';

circH1 = zeros(nH1,1);
validPerim = perimH1 > 0;
circH1(validPerim) = 4*pi*areasH1(validPerim) ./ (perimH1(validPerim).^2);

touchesBorderH1 = false(nH1,1);
for k = 1:nH1
    bb = statsH1(k).BoundingBox;
    x1 = bb(1);
    y1 = bb(2);
    x2 = bb(1) + bb(3);
    y2 = bb(2) + bb(4);

    touchesBorderH1(k) = ...
        x1 <= borderMargin || ...
        y1 <= borderMargin || ...
        x2 >= imgW - borderMargin + 1 || ...
        y2 >= imgH - borderMargin + 1;
end

isNucleus = ...
    areasH1 >= minNucleusArea & ...
    areasH1 <= maxNucleusArea & ...
    solH1 >= minSolidityNucleus & ...
    circH1 >= minCircularityNucleus & ...
    ~touchesBorderH1;

nucleusIdx = find(isNucleus);
nNuclei = numel(nucleusIdx);

if nNuclei == 0
    error('No nuclei detected. Try lowering minNucleusArea or adaptiveSensitivityH1.');
end

nucleusCentroids = reshape([statsH1(nucleusIdx).Centroid], 2, [])';

%% ---------------- BUILD CANDIDATE MN LIST ----------------
candidateCentroids = [];
candidateAreas = [];
candidateCircularity = [];
candidateSolidity = [];
candidateHasEmerin = [];
candidateSource = [];      % 1 = H1 candidate, 2 = bright round Emerin candidate
candidatePixelIdxList = {}; % store pixels so we can debug / deduplicate if needed

%% 1) Standard H1-based MN candidates
nMN = numel(statsMN);

for k = 1:nMN
    A = statsMN(k).Area;
    S = statsMN(k).Solidity;
    P = statsMN(k).Perimeter;

    if P > 0
        C = 4*pi*A/(P^2);
    else
        C = 0;
    end

    bb = statsMN(k).BoundingBox;
    x1 = bb(1); y1 = bb(2);
    x2 = bb(1) + bb(3); y2 = bb(2) + bb(4);

    touchesBorder = ...
        x1 <= borderMargin || ...
        y1 <= borderMargin || ...
        x2 >= imgW - borderMargin + 1 || ...
        y2 >= imgH - borderMargin + 1;

    pix = statsMN(k).PixelIdxList;
    hasEmerin = any(BW_E_support(pix));

    if A >= minMicronucleusArea && A <= maxMicronucleusArea && ...
       S >= minSolidityMN && C >= minCircularityMN && ...
       ~touchesBorder

        candidateCentroids(end+1, :) = statsMN(k).Centroid; %#ok<AGROW>
        candidateAreas(end+1, 1) = A; %#ok<AGROW>
        candidateCircularity(end+1, 1) = C; %#ok<AGROW>
        candidateSolidity(end+1, 1) = S; %#ok<AGROW>
        candidateHasEmerin(end+1, 1) = hasEmerin; %#ok<AGROW>
        candidateSource(end+1, 1) = 1; %#ok<AGROW>
        candidatePixelIdxList{end+1,1} = pix; %#ok<AGROW>
    end
end

%% 2) Extra bright, round Emerin-circle candidates
if useBrightRoundEmerinCandidates
    nE = numel(statsE);

    for k = 1:nE
        A = statsE(k).Area;
        S = statsE(k).Solidity;
        P = statsE(k).Perimeter;
        M = statsE(k).MeanIntensity;

        if P > 0
            C = 4*pi*A/(P^2);
        else
            C = 0;
        end

        bb = statsE(k).BoundingBox;
        x1 = bb(1); y1 = bb(2);
        x2 = bb(1) + bb(3); y2 = bb(2) + bb(4);

        touchesBorder = ...
            x1 <= borderMargin || ...
            y1 <= borderMargin || ...
            x2 >= imgW - borderMargin + 1 || ...
            y2 >= imgH - borderMargin + 1;

        pix = statsE(k).PixelIdxList;
        hasAnyH1 = any(BW_H1_cand(pix));

        if A >= minBrightEmerinArea && A <= maxBrightEmerinArea && ...
           S >= minBrightEmerinSolidity && ...
           C >= minBrightEmerinCircularity && ...
           M >= minBrightEmerinMeanIntensity && ...
           hasAnyH1 && ~touchesBorder

            candidateCentroids(end+1, :) = statsE(k).Centroid; %#ok<AGROW>
            candidateAreas(end+1, 1) = A; %#ok<AGROW>
            candidateCircularity(end+1, 1) = C; %#ok<AGROW>
            candidateSolidity(end+1, 1) = S; %#ok<AGROW>
            candidateHasEmerin(end+1, 1) = true; %#ok<AGROW>
            candidateSource(end+1, 1) = 2; %#ok<AGROW>
            candidatePixelIdxList{end+1,1} = pix; %#ok<AGROW>
        end
    end
end

%% ---------------- REMOVE DUPLICATE CANDIDATES ----------------
% Keep brighter/more supported candidates by preferring:
% 1) Emerin-supported over not
% 2) larger area if near-identical positions
nCandidatesRaw = size(candidateCentroids, 1);
keepCandidate = true(nCandidatesRaw,1);

for i = 1:nCandidatesRaw
    if ~keepCandidate(i)
        continue;
    end

    for j = i+1:nCandidatesRaw
        if ~keepCandidate(j)
            continue;
        end

        d = norm(candidateCentroids(i,:) - candidateCentroids(j,:));

        if d < candidateMergeDistance
            score_i = 10*candidateHasEmerin(i) + candidateAreas(i);
            score_j = 10*candidateHasEmerin(j) + candidateAreas(j);

            if score_i >= score_j
                keepCandidate(j) = false;
            else
                keepCandidate(i) = false;
                break;
            end
        end
    end
end

candidateCentroids = candidateCentroids(keepCandidate, :);
candidateAreas = candidateAreas(keepCandidate);
candidateCircularity = candidateCircularity(keepCandidate);
candidateSolidity = candidateSolidity(keepCandidate);
candidateHasEmerin = candidateHasEmerin(keepCandidate);
candidateSource = candidateSource(keepCandidate);
candidatePixelIdxList = candidatePixelIdxList(keepCandidate);
%% ---------------- MEASURE EMERIN INTENSITY FOR EACH CANDIDATE ----------------
candidateEmerinMeanIntensity = zeros(size(candidateAreas));

for k = 1:numel(candidatePixelIdxList)
    pix = candidatePixelIdxList{k};
    candidateEmerinMeanIntensity(k) = mean(Emerin_enh(pix));
end

%% ---------------- ASSIGN CANDIDATES TO NEAREST NUCLEUS ----------------
nCandidates = size(candidateCentroids, 1);

candidateMN = true(nCandidates,1);
isMN = false(nCandidates,1);
assignedParent = nan(nCandidates,1);
distanceToParent = nan(nCandidates,1);
insideAnyNucleus = false(nCandidates,1);
rescuedBrightIntranuclearMN = false(nCandidates,1);

% Build one combined mask of all detected nuclei
allNucleiMask = false(imgH, imgW);
for ii = 1:numel(nucleusIdx)
    objIdx = nucleusIdx(ii);
    allNucleiMask(statsH1(objIdx).PixelIdxList) = true;
end

for k = 1:nCandidates
    cMN = candidateCentroids(k, :);

    % nearest nucleus assignment still used for parent relationship
    dists = sqrt(sum((nucleusCentroids - cMN).^2, 2));
    [dmin, p] = min(dists);

    if dmin <= maxDistanceToNearestNucleus
        assignedParent(k) = p;
        distanceToParent(k) = dmin;

        % Check whether MN centroid lies inside ANY nucleus
        mn_x = round(cMN(1));
        mn_y = round(cMN(2));

        mn_x = max(1, min(imgW, mn_x));
        mn_y = max(1, min(imgH, mn_y));

        if allNucleiMask(mn_y, mn_x)
            insideAnyNucleus(k) = true;
        end

        % Default behavior: keep if not inside nucleus
        if ~insideAnyNucleus(k)
            isMN(k) = true;

        else
            % Rescue rule: very bright, round, compact emerin object
            if allowVeryBrightIntranuclearEmerinMN
                if candidateEmerinMeanIntensity(k) >= minRescueEmerinMeanIntensity && ...
                   candidateCircularity(k) >= minRescueEmerinCircularity && ...
                   candidateSolidity(k) >= minRescueEmerinSolidity && ...
                   candidateAreas(k) >= minRescueEmerinArea && ...
                   candidateAreas(k) <= maxRescueEmerinArea

                    isMN(k) = true;
                    rescuedBrightIntranuclearMN(k) = true;
                end
            end
        end
    end
end

%% ---------------- REMOVE NUCLEI WITH >=3 MN ----------------
mnCounts = zeros(nNuclei,1);

for p = 1:nNuclei
    mnCounts(p) = sum((assignedParent == p) & isMN);
end

badParents = find(mnCounts >= 3);

for i = 1:numel(badParents)
    p = badParents(i);
    isMN((assignedParent == p) & isMN) = false;
end

goodNuclei = true(nNuclei,1);
goodNuclei(badParents) = false;

mnCounts = zeros(nNuclei,1);
for p = 1:nNuclei
    if goodNuclei(p)
        mnCounts(p) = sum((assignedParent == p) & isMN);
    else
        mnCounts(p) = 0;
    end
end

%% ---------------- FINAL COUNTS ----------------
nValidNuclei = sum(goodNuclei);
nMicronuclei = sum(isMN);
nMicronucleatedNuclei = numel(unique(assignedParent(isMN)));

nMN_H1_only = sum(isMN & ~candidateHasEmerin);
nMN_H1_and_Emerin = sum(isMN & candidateHasEmerin);
nMN_brightRoundEmerinSource = sum(isMN & candidateSource == 2);

fprintf('Originally detected nuclei: %d\n', nNuclei);
fprintf('Valid nuclei after >=3 MN exclusion: %d\n', nValidNuclei);
fprintf('Micronuclei total: %d\n', nMicronuclei);
fprintf('Micronucleated nuclei: %d\n', nMicronucleatedNuclei);
fprintf('MN that are H1 only: %d\n', nMN_H1_only);
fprintf('MN that are H1 + Emerin: %d\n', nMN_H1_and_Emerin);
fprintf('MN added from bright round Emerin candidates: %d\n', nMN_brightRoundEmerinSource);

%% ---------------- SAVE TABLES ----------------
nucleusSummary = table( ...
    (1:nNuclei)', ...
    nucleusIdx(:), ...
    nucleusCentroids(:,1), ...
    nucleusCentroids(:,2), ...
    mnCounts, ...
    ~goodNuclei, ...
    'VariableNames', ...
    {'NucleusLabel','H1_ObjectIndex','CentroidX','CentroidY','NumAssignedMN','ExcludedBecause3OrMoreMN'});

writetable(nucleusSummary, 'nucleus_mn_summary.csv');

mnIdx = find(isMN);

if ~isempty(mnIdx)
    mnSummary = table( ...
        mnIdx(:), ...
        candidateCentroids(mnIdx,1), ...
        candidateCentroids(mnIdx,2), ...
        candidateAreas(mnIdx), ...
        candidateCircularity(mnIdx), ...
        candidateSolidity(mnIdx), ...
        candidateEmerinMeanIntensity(mnIdx), ...
        candidateHasEmerin(mnIdx), ...
        candidateSource(mnIdx), ...
        assignedParent(mnIdx), ...
        distanceToParent(mnIdx), ...
        insideAnyNucleus(mnIdx), ...
        rescuedBrightIntranuclearMN(mnIdx), ...
        'VariableNames', ...
        {'MN_ObjectIndex','CentroidX','CentroidY','Area','Circularity','Solidity', ...
         'EmerinMeanIntensity','HasEmerinSupport','Source_1H1_2BrightRoundEmerin', ...
         'AssignedParentNucleus','DistanceToParent','InsideAnyNucleus', ...
         'RescuedBrightIntranuclearMN'});

    writetable(mnSummary, 'micronucleus_parent_assignments.csv');
end
%% ---------------- DISPLAY FINAL OVERLAY ON EMERIN ----------------
figure;
imshow(Emerin, []); hold on;
title('Green = valid nuclei, Red = final MN on Emerin image');

for i = find(goodNuclei)'
    c = statsH1(nucleusIdx(i)).Centroid;
    plot(c(1), c(2), 'go', 'MarkerSize', 4, 'LineWidth', 1);
end

for k = find(isMN)'
    c = candidateCentroids(k,:);
    plot(c(1), c(2), 'ro', 'MarkerSize', 4, 'LineWidth', 1);
end

hold off;

%% ---------------- DISPLAY FINAL OVERLAY ON H1 ----------------
figure;
imshow(H1, []); hold on;
title('Green = valid nuclei, Red = final MN on H1 image');

% Plot nuclei
for i = find(goodNuclei)'
    c = statsH1(nucleusIdx(i)).Centroid;
    plot(c(1), c(2), 'go', 'MarkerSize', 4, 'LineWidth', 1);
end

% Plot MN
for k = find(isMN)'
    c = candidateCentroids(k,:);
    plot(c(1), c(2), 'ro', 'MarkerSize', 4, 'LineWidth', 1);
end

hold off;

%% ---------------- OPTIONAL DEBUG FIGURE ON EMERIN ----------------
if showDebugFigures
    figure;
    imshow(Emerin, []); hold on;
    title('Yellow = candidate MN, Red = final MN, Green = valid nuclei');

    for k = find(candidateMN)'
        c = candidateCentroids(k,:);
        plot(c(1), c(2), 'yo', 'MarkerSize', 4, 'LineWidth', 1);
    end

    for k = find(isMN)'
        c = candidateCentroids(k,:);
        plot(c(1), c(2), 'ro', 'MarkerSize', 5, 'LineWidth', 1.2);
    end

    for i = find(goodNuclei)'
        c = statsH1(nucleusIdx(i)).Centroid;
        plot(c(1), c(2), 'go', 'MarkerSize', 5, 'LineWidth', 1.2);
    end

    hold off;
end