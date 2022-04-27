function varargout = GUI2(varargin)
% GUI2 MATLAB code for GUI2.fig
%      GUI2, by itself, creates a new GUI2 or raises the existing
%      singleton*.
%
%      H = GUI2 returns the handle to a new GUI2 or the handle to
%      the existing singleton*.
%
%      GUI2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI2.M with the given input arguments.
%
%      GUI2('Property','Value',...) creates a new GUI2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI2

% Last Modified by GUIDE v2.5 28-Mar-2022 11:18:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI2_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI2_OutputFcn, ...
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


% --- Executes just before GUI2 is made visible.
function GUI2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI2 (see VARARGIN)

% Choose default command line output for GUI2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI2_OutputFcn(hObject, eventdata, handles) 
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
f1 = 30;
fai1 = pi/3;
Fs = 2^9;
T =1;      % T = N/Fs 
t = 0:1/Fs:T-1/Fs;
N = length(t);
NFFT = 2^nextpow2(N);%nextpow2(N),靠的最近的2的指数
x = sin(2*pi*f1*t);
% ======================= fft ==================================
deltaF = 1/T;   %  deltaF = dFs/(N-1);%如果是表格或mat文件，推荐用这一种
vecf = (0:N-1)*deltaF;
% vecf  = linspace(0,Fs,N); 
% linspace(x1,x2,n) % 生成 n 个点。这些点的间距为 (x2-x1)/(n-1)。
% 如果数据x是奇数，f=vecf
tic;
xk = 2*fft(x,NFFT)/N;        
toc;

Ampli = abs(xk);
phase = unwrap(angle(xk));
subplot(3,1,1);plot(t,x);title('信号')
subplot(3,1,2);plot(vecf,Ampli);title('FFT双边谱')
subplot(3,1,3);plot(vecf(1:Fs/2+1),Ampli(1:Fs/2+1));title('FFT单边谱')



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 功能：离散时间傅里叶变换DTFT
N=8;
%原离散信号有8点
n=[0:1:N-1];
%原信号是1行8列的矩阵
xn=0.5.^n;
%构建原始信号，为指数信号
% ========== 求dtft变换，采用原始定义对复指数分量求和 ============
w=[-800:1:800]*4*pi/800;
%频域共-800--+800的长度
%本应是无穷，高频分量很少，故省去
martrix = w'.*n;
X = exp(-1i*(martrix))*xn';
% ======================= figure ==========================
subplot(211)
stem(n,xn);
title('原始信号(指数信号)');
subplot(212);
stem(w/pi,abs(X));
title('DTFT变换')


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;
clc;
close GUI2;