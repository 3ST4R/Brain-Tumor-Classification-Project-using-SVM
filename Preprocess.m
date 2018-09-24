function [ resultimg ] = Preprocess( img )
[M N] = size(img);
n = zeros(M+2,N+2);
med = zeros(M+2,N+2);
[R C] = size(n);
n(2:R-1,2:C-1) = img(1:end,1:end);
med(2:R-1,2:C-1) = img(1:end,1:end);
for i = 2:R-1
    for j = 2:C-1
        temp = [med(i-1,j-1) med(i-1,j) med(i-1,j+1);med(i,j-1) med(i,j) med(i,j+1);....
                med(i+1,j-1) med(i+1,j) med(i+1,j+1)];
            medsort = sort(temp(:),'ascend');
            med(i,j) = medsort(5);
    end
end
demed = med(2:R-1,2:C-1);
resultimg = uint8(demed);

end

