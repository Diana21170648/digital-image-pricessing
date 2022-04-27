function varargout = GUI4(varargin)
% GUI4 MATLAB code for GUI4.fig
%      GUI4, by itself, creates a new GUI4 or raises the existing
%      singleton*.
%
%      H = GUI4 returns the handle to a new GUI4 or the handle to
%      the existing singleton*.
%
%      GUI4('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI4.M with the given input arguments.
%
%      GUI4('Property','Value',...) creates a new GUI4 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI4_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI4_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI4

% Last Modified by GUIDE v2.5 05-Apr-2022 22:10:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI4_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI4_OutputFcn, ...
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


% --- Executes just before GUI4 is made visible.
function GUI4_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI4 (see VARARGIN)

% Choose default command line output for GUI4
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI4 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI4_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%P-Hough --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I=imread('车道检测.jpg');     %读图
rotI=rgb2gray(I);               %rgb转灰度图
subplot(2,2,1);               %显示灰度图
imshow(rotI);
title('灰度图像');
axis on;
BW=edge(rotI,'prewitt');
subplot(2,2,2);
imshow(BW);
title('prewitt算子边缘检测 后图像');
axis on;
[H,T,R]=hough(BW);
subplot(2,2,3);
imshow(H,[],'XData',T,'YData',R,'InitialMagnification','fit');
title('霍夫变换图');
xlabel('\theta'),ylabel('\rho');
axis on , axis normal, hold on;
P=houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
x=T(P(:,2));y=R(P(:,1));
plot(x,y,'s','color','white');
lines=houghlines(BW,T,R,P,'FillGap',5,'MinLength',7);
subplot(2,2,4);,imshow(rotI);
title('霍夫变换图像检测');
axis on;
hold on;
max_len=0;
for k=1:length(lines)
xy=[lines(k).point1;lines(k).point2];
plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
len=norm(lines(k).point1-lines(k).point2);
if(len>max_len)
max_len=len;
xy_long=xy;
end
end
plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','cyan');
