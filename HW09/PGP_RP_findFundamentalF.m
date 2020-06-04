%%original:by Stepan Pazekha
%%optimized: PGP
function F = PGP_RP_findFundamentalF(PIleftfname, PIrightfname, PEpileftfname, PEpirightfname)
% function F = PGP_findFundamentalF(I_left, I_right, P_1, P_2)
% PIleftfname, PIrightfname - stereo image pair.
% Pcorrfnamr: file name of corresponding points from stereo images,
%  must have dimensions of [n, 2].

%read images
F = 0;
[h, w, c] = size(imread(PIleftfname));
if(c == 3) I_left = rgb2gray(imread(PIleftfname)); end;

[h2, w2, c] = size(imread(PIrightfname));
if(c == 3) I_right = rgb2gray(imread(PIrightfname)); end;

w = min(w, w2);
h = min(h, h2);
I_left = I_left(1:h, 1:w);
I_right = I_right(1:h, 1:w);

%read correspondances
epipolesLEFTpts=importdata(PEpileftfname);
epipolesRIGHTpts=importdata(PEpirightfname);
%very small sanity check
[nL, cL] = size(epipolesLEFTpts);
[nR, cR] = size(epipolesRIGHTpts);
if (nL ~= nR)
    error('Error! number of correspondances have to match');
    return;
end;

% CALCULATE F
F = PGP_RP_calcFundamentalF(epipolesLEFTpts, epipolesRIGHTpts);
F = F/F(3,3);
%
disp('Fundamental matrix: ');
disp(F);
[U SS V] = svd(F);
el = V(:,3) / V(3,3);
er = U(:,3) / U(3,3);
disp('Epipole in the left image: ');
disp(el);
disp('Epipole in the right image: ');
disp(er);
%
%calculate Epipol lines
%epiLEFTlines = transpose( F * [transpose(epipolesLEFTpts);ones(1,nL)]);
%epiRIGHTlines = transpose(transpose(F) * [transpose(epipolesRIGHTpts);ones(1,nR)]);
epiLEFTlines = [epipolesLEFTpts ones(nL,1)] * transpose(F)
epiRIGHTlines = [epipolesRIGHTpts ones(nR,1)] * F
%calculated begin / end coordinates(in image) for those epipol lines
%first image
leftXs = ones(1,nL);
rightXs = w*ones(1,nL);
leftYs=PGPInHom2Ycoord(epiLEFTlines,1);
rightYs=PGPInHom2Ycoord(epiLEFTlines,w);
%draw
colormap gray;
subplot(1, 2, 1);
imagesc(I_left); hold on;
scatter(epipolesLEFTpts(:,1),epipolesLEFTpts(:,2),1000,'rx','LineWidth',3)
allX=[leftXs;rightXs];
allY=[leftYs';rightYs'];
line(allX,allY, 'Color', [0, 1, 0]); hold off;
% second image, draw
subplot(1, 2, 2);
leftYs=PGPInHom2Ycoord(epiRIGHTlines,1);
rightYs=PGPInHom2Ycoord(epiRIGHTlines,w);
imagesc(I_right); hold on;
scatter(epipolesRIGHTpts(:,1),epipolesRIGHTpts(:,2),1000,'rx','LineWidth',3)
allX=[leftXs;rightXs];
allY=[leftYs';rightYs'];
line(allX,allY, 'Color', [0, 1, 0]); hold off;
end
