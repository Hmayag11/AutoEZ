
rawdata = load('ct.mat');
J = permute(rawdata.ct.Cube, [1 2 4 3]);
dicomwrite(J,'yourFile.dcm');
dicomanon('yourFile.dcm','CT_P31.dcm')
X = dicomread('CT_P31.dcm');
 info = dicominfo('CT_P31.dcm')
%%
 addpath('C:\Users\hpartamian-temp\Downloads\xiangruili-dicm2nii-3fe1a27')
 dicm2nii

 %% Load DICOM 
 filename='MRI_P1.dcm';
 X = dicomread(filename);
 % PLOT A SLICE,
 s=150;%slice number
 figure;imagesc(X(:,:,1,150))
 info = dicominfo( filename);
 info.PatientName
 % plot in 3D
 sx = 1;
 sy = 1;
 sz = 2.5;
 A = [sx 0 0 0; 0 sy 0 0; 0 0 sz 0; 0 0 0 1];
 tform = affinetform3d(A);
 J = permute(X, [1 2 4 3]);
 vol = volshow(squeeze(J(:,:,:,1)),renderingStyle="SlicePlanes",Transformation=tform);