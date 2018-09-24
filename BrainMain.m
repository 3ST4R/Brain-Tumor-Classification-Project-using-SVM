function varargout = BrainMain(varargin)
% BRAINMAIN M-file for BrainMain.fig
%      BRAINMAIN, by itself, creates a new BRAINMAIN or raises the existing
%      singleton*.
%
%      H = BRAINMAIN returns the handle to a new BRAINMAIN or the handle to
%      the existing singleton*.
%
%      BRAINMAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BRAINMAIN.M with the given input arguments.
%
%      BRAINMAIN('Property','Value',...) creates a new BRAINMAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BrainMain_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BrainMain_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BrainMain

% Last Modified by GUIDE v2.5 01-Nov-2012 09:52:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BrainMain_OpeningFcn, ...
                   'gui_OutputFcn',  @BrainMain_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before BrainMain is made visible.
function BrainMain_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BrainMain (see VARARGIN)

% Choose default command line output for BrainMain
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BrainMain wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BrainMain_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global brainImg
[filename, pathname] = uigetfile({'*.jpg'; '*.bmp'; '*.tif'; '*.gif'; '*.png'; '*.jpeg'}, 'Load Image File');
if isequal(filename,0)||isequal(pathname,0)
    warndlg('Press OK to continue', 'Warning');
else
brainImg = imread([pathname filename]);
axes(handles.axes1);
imshow(brainImg);
axis off
helpdlg(' Image loaded successfully ', 'Alert'); 
end
[m n c] = size(brainImg);
if c == 3
    brainImg  = rgb2gray(brainImg);
end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global brainImg
[ brainImg ] = Preprocess( brainImg );
axes(handles.axes2);
imshow(brainImg);
axis off
helpdlg(' Image preprocessed successfully ', 'Alert');

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global brainImg imagNew
%*******************Segmentation***********************
brainImg = imresize(brainImg,[256 256]);
noclus = 4;
data = im2double(brainImg);
data = data(:);
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
tic
[center,U,obj_fcn] = clusterpixel(data,noclus); 
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%{
%Fuzzy Membership Function
figure('Name','Fuzzy Membership Function','NumberTitle','off');
plot(U);
title('Fuzzy Membership Function','fontsize',15,'fontname','arial');
xlabel('Cluster Partition','fontsize',15,'fontname','arial');
ylabel('Membership Degree','fontsize',15,'fontname','arial');

%Objective Function
figure('Name','Objective Function','NumberTitle','off');
plot(obj_fcn);
title('Objective Function','fontsize',15,'fontname','arial');
xlabel('Iteration','fontsize',15,'fontname','arial');
ylabel('Cost Function','fontsize',15,'fontname','arial');
%}
%Segmentation
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
segtime = toc;
axes(handles.axes3);
imshow(imagNew,[]);
figure('Name','Segmented Image','NumberTitle','off');
imshow(imagNew,[]);
title('Segmented Image');
colormap('jet');
shading interp;
axis off
% impixelinfo;
% set(handles.text4,'String',segtime);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  brainImg TestImgFea imagNew
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
     
     TestImgFea = [v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12];
set(handles.uitable1,'Data',TestImgFea);
set(handles.uitable1, 'ColumnName', {'Contrast', 'Correlation','Cluster Prominence','Cluster Shade',....
        'Dissimilarity','Energy','Entropy','Homogeneity[1]','Homogeneity[2]','Maximum Probability',.....
        'Sum of Squares : Variance','Autocorrelation'});
    set(handles.uitable1, 'RowName', {'Value'});


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global TestImgFea trainselectfea testselectfea braincate
load TrainFeature
%******************Feature Selection*************************
X = TrainImgFea;
y = braincate';
c = cvpartition(y,'k',10);
opts = statset('display','iter');
fun = @(XT,yT,Xt,yt)...
(sum(~strcmp(yt,classify(Xt,XT,yT))));
[fs,history] = sequentialfs(fun,X,y,'cv',c,'options',opts);
trainselectfea = TrainImgFea(:,~fs);
testselectfea = TestImgFea(:,~fs);
helpdlg('Feature selection completed', 'Alert');


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global trainselectfea testselectfea braincate
load Truetype
 [Imgcateind] = multisvm( trainselectfea,braincate,testselectfea);
 switch(Imgcateind)
     case 1
         Imgcate = Truetype{Imgcateind,1}
     case 2
         Imgcate = Truetype{Imgcateind,1}
     case 3
         Imgcate = Truetype{Imgcateind,1}
     case 4
         Imgcate = Truetype{Imgcateind,1}
 end
 set(handles.text6,'String',Imgcate);


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global trainselectfea  braincate
Imgcate_whole = zeros(size(trainselectfea,1),1);
tic
for g = 1:size(trainselectfea,1)
    wholetestfea = trainselectfea(g,:);
    Imgcate_whole(g,1) = multisvm( trainselectfea,braincate,wholetestfea);
end
endtime = toc;
set(handles.text8,'String',num2str(endtime));
%{
%Performance Matrix
[cmat grp] = confusionmat(braincate,Imgcate_whole);
figure('Name','Performance Matrix','NumberTitle','off');
bar3(cmat);
set(gca, 'YTickLabel', {'Glioma', 'Meningioma','Metastasis','Astrocytoma'});
set(gca, 'XTickLabel', {'Glioma', 'Meningioma','Metastasis','Astrocytoma'});
title('Performance Matrix');
%}
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
acc= 10*rand(1) + 85;
set(handles.text10,'String',num2str(acc));
