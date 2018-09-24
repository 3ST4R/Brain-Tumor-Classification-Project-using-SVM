TrainImgFea = zeros(39,12);
for i = 1:39
    trainimg = imread(['Images\',num2str(i),'.jpg']);
    imshow(trainimg);
    title(['Image- ',num2str(i)]);
    [r c d] = size(trainimg);
    if d == 3
        trainimg = rgb2gray(trainimg);
    end
    trainimg = Preprocess(trainimg);
    %^^^^^^^^^^^^^^Segmentation^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    %*******************Segmentation***********************
    trainimg = imresize(trainimg,[256 256]);
    noclus = 4;
    data = im2double(trainimg);
    data = data(:);
    %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      [center,U,obj_fcn] = clusterpixel(data,noclus); 
    %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    maxU = max(U);
    index1 = find(U(1,:) == maxU);
    index2 = find(U(2,:) == maxU);
    index3 = find(U(3,:) == maxU);
    index4 = find(U(4,:) == maxU);
    fcmImage(1:length(data))=0;       
    fcmImage(index1)= 1;
    fcmImage(index2)= 0.8;
    fcmImage(index3)= 0.6;
    fcmImage(index4)= 0.4;
    imagNew = reshape(fcmImage,256,256);
    
    %^^^^^^^^^^^^^^^^^^^^Feature Extraction^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    GLCM_mat = graycomatrix(imagNew,'Offset',[2 0;0 2]);
     
     GLCMstruct = Computefea(GLCM_mat,0);
     
     v1=GLCMstruct.contr(1);

     v2=GLCMstruct.corrm(1);

     v3=GLCMstruct.cprom(1);

     v4=GLCMstruct.cshad(1);

     v5=GLCMstruct.dissi(1);

     v6=GLCMstruct.energ(1);

     v7=GLCMstruct.entro(1);

     v8=GLCMstruct.homom1(1);

     v9=GLCMstruct.homop(1);

     v10=GLCMstruct.maxpr(1);

     v11=GLCMstruct.sosvh(1);

     v12=GLCMstruct.autoc(1);
     
     TrainImgFea(i,:) = [v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12];
    %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     
end
braincate(1:8) = 1;
braincate(9:16) = 2;
braincate(17:34) = 3;
braincate(35:39) = 4;
save TrainFeature 'TrainImgFea' 'braincate'
Truetype{1,1} = 'Glioma';
ruetype{2,1} = 'Meningioma';
Truetype{3,1} = 'Metastasis';
Truetype{4,1} = 'Astrocytoma';
% for j = 1:8
%     Truetype{j,1} = 'Glioma';
% end
% for j = 9:16
%     Truetype{j,1} = 'Meningioma';
% end
% for j = 17:34
%     Truetype{j,1} = 'Metastasis';
% end
% for j = 35:39
%     Truetype{j,1} = 'Astrocytoma';
% end
save Truetype Truetype