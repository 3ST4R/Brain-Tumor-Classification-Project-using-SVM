function [feastruct] = Computefea(glcmin,pairs)

% If 'pairs' not entered: set pairs to 0 
if ((nargin > 2) || (nargin == 0))
   error('Too many or too few input arguments. Enter GLCM and pairs.');
elseif ( (nargin == 2) ) 
    if ((size(glcmin,1) <= 1) || (size(glcmin,2) <= 1))
       error('The GLCM should be a 2-D or 3-D matrix.');
    elseif ( size(glcmin,1) ~= size(glcmin,2) )
        error('Each GLCM should be square with NumLevels rows and NumLevels cols');
    end    
elseif (nargin == 1) % only GLCM is entered
    pairs = 0; % default is numbers and input 1 for percentage
    if ((size(glcmin,1) <= 1) || (size(glcmin,2) <= 1))
       error('The GLCM should be a 2-D or 3-D matrix.');
    elseif ( size(glcmin,1) ~= size(glcmin,2) )
       error('Each GLCM should be square with NumLevels rows and NumLevels cols');
    end    
end


format long e
if (pairs == 1)
    newn = 1;
    for nglcm = 1:2:size(glcmin,3)
        glcm(:,:,newn)  = glcmin(:,:,nglcm) + glcmin(:,:,nglcm+1);
        newn = newn + 1;
    end
elseif (pairs == 0)
    glcm = glcmin;
end

size_glcm_1 = size(glcm,1);
size_glcm_2 = size(glcm,2);
size_glcm_3 = size(glcm,3);

% checked 
feastruct.autoc = zeros(1,size_glcm_3); % Autocorrelation 
feastruct.contr = zeros(1,size_glcm_3); % Contrast
feastruct.corrm = zeros(1,size_glcm_3); % Correlation
feastruct.corrp = zeros(1,size_glcm_3); % Correlation
feastruct.cprom = zeros(1,size_glcm_3); % Cluster Prominence
feastruct.cshad = zeros(1,size_glcm_3); % Cluster Shade
feastruct.dissi = zeros(1,size_glcm_3); % Dissimilarity
feastruct.energ = zeros(1,size_glcm_3); % Energy
feastruct.entro = zeros(1,size_glcm_3); % Entropy
feastruct.homom1 = zeros(1,size_glcm_3); % Homogeneity
feastruct.homop = zeros(1,size_glcm_3); % Homogeneity
feastruct.maxpr = zeros(1,size_glcm_3); % Maximum probability

feastruct.sosvh = zeros(1,size_glcm_3); % Sum of sqaures
feastruct.savgh = zeros(1,size_glcm_3); % Sum average 
feastruct.svarh = zeros(1,size_glcm_3); % Sum variance 
feastruct.senth = zeros(1,size_glcm_3); % Sum entropy 
feastruct.dvarh = zeros(1,size_glcm_3); % Difference variance 
%feastruct.dvarh2 = zeros(1,size_glcm_3); % Difference variance 
feastruct.denth = zeros(1,size_glcm_3); % Difference entropy 
feastruct.inf1h = zeros(1,size_glcm_3); % Information measure of correlation1 
feastruct.inf2h = zeros(1,size_glcm_3); % Informaiton measure of correlation2 
%feastruct.mxcch = zeros(1,size_glcm_3);% maximal correlation coefficient 
%feastruct.invdc = zeros(1,size_glcm_3);% Inverse difference (INV)  
feastruct.indnc = zeros(1,size_glcm_3); % Inverse difference normalized 
feastruct.idmnc = zeros(1,size_glcm_3); % Inverse difference moment normalized 



glcm_sum  = zeros(size_glcm_3,1);
glcm_mean = zeros(size_glcm_3,1);
glcm_var  = zeros(size_glcm_3,1);


u_x = zeros(size_glcm_3,1);
u_y = zeros(size_glcm_3,1);
s_x = zeros(size_glcm_3,1);
s_y = zeros(size_glcm_3,1);




p_x = zeros(size_glcm_1,size_glcm_3); 
p_y = zeros(size_glcm_2,size_glcm_3); 
p_xplusy = zeros((size_glcm_1*2 - 1),size_glcm_3); 
p_xminusy = zeros((size_glcm_1),size_glcm_3);

hxy  = zeros(size_glcm_3,1);
hxy1 = zeros(size_glcm_3,1);
hx   = zeros(size_glcm_3,1);
hy   = zeros(size_glcm_3,1);
hxy2 = zeros(size_glcm_3,1);



for k = 1:size_glcm_3 % number glcms

    glcm_sum(k) = sum(sum(glcm(:,:,k)));
    glcm(:,:,k) = glcm(:,:,k)./glcm_sum(k); % Normalize each glcm
    glcm_mean(k) = mean2(glcm(:,:,k)); % compute mean after norm
    glcm_var(k)  = (std2(glcm(:,:,k)))^2;
    
    for i = 1:size_glcm_1

        for j = 1:size_glcm_2

            feastruct.contr(k) = feastruct.contr(k) + (abs(i - j))^2.*glcm(i,j,k);
            feastruct.dissi(k) = feastruct.dissi(k) + (abs(i - j)*glcm(i,j,k));
            feastruct.energ(k) = feastruct.energ(k) + (glcm(i,j,k).^2);
            feastruct.entro(k) = feastruct.entro(k) - (glcm(i,j,k)*log(glcm(i,j,k) + eps));
            feastruct.homom1(k) = feastruct.homom1(k) + (glcm(i,j,k)/( 1 + abs(i-j) ));
            feastruct.homop(k) = feastruct.homop(k) + (glcm(i,j,k)/( 1 + (i - j)^2));
           
            feastruct.sosvh(k) = feastruct.sosvh(k) + glcm(i,j,k)*((i - glcm_mean(k))^2);
            
            %feastruct.invdc(k) = feastruct.homom1(k);
            feastruct.indnc(k) = feastruct.indnc(k) + (glcm(i,j,k)/( 1 + (abs(i-j)/size_glcm_1) ));
            feastruct.idmnc(k) = feastruct.idmnc(k) + (glcm(i,j,k)/( 1 + ((i - j)/size_glcm_1)^2));
            u_x(k)          = u_x(k) + (i)*glcm(i,j,k);
            u_y(k)          = u_y(k) + (j)*glcm(i,j,k);
            
        end
        
    end
    feastruct.maxpr(k) = max(max(glcm(:,:,k)));
end

for k = 1:size_glcm_3
    
    for i = 1:size_glcm_1
        
        for j = 1:size_glcm_2
            p_x(i,k) = p_x(i,k) + glcm(i,j,k);
            p_y(i,k) = p_y(i,k) + glcm(j,i,k); % taking i for j and j for i
            if (ismember((i + j),[2:2*size_glcm_1])) 
                p_xplusy((i+j)-1,k) = p_xplusy((i+j)-1,k) + glcm(i,j,k);
            end
            if (ismember(abs(i-j),[0:(size_glcm_1-1)])) 
                p_xminusy((abs(i-j))+1,k) = p_xminusy((abs(i-j))+1,k) +...
                    glcm(i,j,k);
            end
        end
    end
    

    
end


for k = 1:(size_glcm_3)
    
    for i = 1:(2*(size_glcm_1)-1)
        feastruct.savgh(k) = feastruct.savgh(k) + (i+1)*p_xplusy(i,k);
        % the summation for savgh is for i from 2 to 2*Ng hence (i+1)
        feastruct.senth(k) = feastruct.senth(k) - (p_xplusy(i,k)*log(p_xplusy(i,k) + eps));
    end

end
% compute sum variance with the help of sum entropy
for k = 1:(size_glcm_3)
    
    for i = 1:(2*(size_glcm_1)-1)
        feastruct.svarh(k) = feastruct.svarh(k) + (((i+1) - feastruct.senth(k))^2)*p_xplusy(i,k);
        % the summation for savgh is for i from 2 to 2*Ng hence (i+1)
    end

end
% compute difference variance, difference entropy, 
for k = 1:size_glcm_3

    for i = 0:(size_glcm_1-1)
        feastruct.denth(k) = feastruct.denth(k) - (p_xminusy(i+1,k)*log(p_xminusy(i+1,k) + eps));
        feastruct.dvarh(k) = feastruct.dvarh(k) + (i^2)*p_xminusy(i+1,k);
    end
end

% compute information measure of correlation
for k = 1:size_glcm_3
    hxy(k) = feastruct.entro(k);
    for i = 1:size_glcm_1
        
        for j = 1:size_glcm_2
            hxy1(k) = hxy1(k) - (glcm(i,j,k)*log(p_x(i,k)*p_y(j,k) + eps));
            hxy2(k) = hxy2(k) - (p_x(i,k)*p_y(j,k)*log(p_x(i,k)*p_y(j,k) + eps));
%             for Qind = 1:(size_glcm_1)
%                 Q(i,j,k) = Q(i,j,k) +...
%                     ( glcm(i,Qind,k)*glcm(j,Qind,k) / (p_x(i,k)*p_y(Qind,k)) ); 
%             end
        end
        hx(k) = hx(k) - (p_x(i,k)*log(p_x(i,k) + eps));
        hy(k) = hy(k) - (p_y(i,k)*log(p_y(i,k) + eps));
    end
    feastruct.inf1h(k) = ( hxy(k) - hxy1(k) ) / ( max([hx(k),hy(k)]) );
    feastruct.inf2h(k) = ( 1 - exp( -2*( hxy2(k) - hxy(k) ) ) )^0.5;
%     eig_Q(k,:)   = eig(Q(:,:,k));
%     sort_eig(k,:)= sort(eig_Q(k,:),'descend');
%     feastruct.mxcch(k) = sort_eig(k,2)^0.5;

end

corm = zeros(size_glcm_3,1);
corp = zeros(size_glcm_3,1);

for k = 1:size_glcm_3
    for i = 1:size_glcm_1
        for j = 1:size_glcm_2
            s_x(k)  = s_x(k)  + (((i) - u_x(k))^2)*glcm(i,j,k);
            s_y(k)  = s_y(k)  + (((j) - u_y(k))^2)*glcm(i,j,k);
            corp(k) = corp(k) + ((i)*(j)*glcm(i,j,k));
            corm(k) = corm(k) + (((i) - u_x(k))*((j) - u_y(k))*glcm(i,j,k));
            feastruct.cprom(k) = feastruct.cprom(k) + (((i + j - u_x(k) - u_y(k))^4)*...
                glcm(i,j,k));
            feastruct.cshad(k) = feastruct.cshad(k) + (((i + j - u_x(k) - u_y(k))^3)*...
                glcm(i,j,k));
        end
    end
    
    s_x(k) = s_x(k) ^ 0.5;
    s_y(k) = s_y(k) ^ 0.5;
    feastruct.autoc(k) = corp(k);
    feastruct.corrp(k) = (corp(k) - u_x(k)*u_y(k))/(s_x(k)*s_y(k));
    feastruct.corrm(k) = corm(k) / (s_x(k)*s_y(k));
%     % alternate values of u and s
%     feastruct.corrp2(k) = (corp(k) - u_x2(k)*u_y2(k))/(s_x2(k)*s_y2(k));
%     feastruct.corrm2(k) = corm(k) / (s_x2(k)*s_y2(k));
end


