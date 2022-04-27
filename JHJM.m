function varargout = JHJM(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @JHJM_OpeningFcn, ...
                   'gui_OutputFcn',  @JHJM_OutputFcn, ...
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
function JHJM_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = JHJM_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function radiobutton1_Callback(hObject, eventdata, handles)

function radiobutton2_Callback(hObject, eventdata, handles)

function radiobutton3_Callback(hObject, eventdata, handles)

function radiobutton4_Callback(hObject, eventdata, handles)

function pushbutton3_Callback(hObject, eventdata, handles)
[filename,pathname,filterindex]=...
uigetfile({'*.*';'*.bmp';'*.tif';'*.png';'*.jpg';'*.jpeg'},'select picture');
str=[pathname filename];  
s=str;
handles.filebig=filterindex;
if filterindex==0
%msgbox('选择图像失败！','error');
return
else   
im=imread(str);  
end 
axes(handles.axes1);  
imshow(im);  %显示图片
handles.img=im;
guidata(hObject,handles);

function popupmenu1_Callback(hObject, eventdata, handles)

function popupmenu1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function radiobutton5_Callback(hObject, eventdata, handles)

 global T
if handles.filebig==0 
%msgbox('请输入函数图像！','error');
return;
else
 axes(handles.axes2);
 imshow(handles.img);
T=handles.img;
axes(handles.axes2);
prompt={'请输入模板维度：'};
defans={'3'};
p=inputdlg(prompt,'input',1,defans);
p1=str2num(p{1});
h1=fspecial('average',[p1 p1]);
I=imfilter(handles.img,h1);
end
imshow(I);
handles.img=I;
guidata(hObject,handles);

 %锐化sobel滤波
function radiobutton6_Callback(hObject, eventdata, handles)
global T
if handles.filebig==0
%msgbox('请输入函数图像！','error');
return;
else
 axes(handles.axes2);
 imshow(handles.img);
T=handles.img;
axes(handles.axes2);
h=fspecial('sobel');
g2=imfilter(handles.img,h);
g3=imadd(g2,handles.img);
end
imshow(g3);
handles.img=g3;
guidata(hObject,handles);

 %锐化prewitt滤波
function radiobutton7_Callback(hObject, eventdata, handles)
global T
if handles.filebig==0
%msgbox('请输入函数图像！','error');
return;
else
 axes(handles.axes2);
 imshow(handles.img);
T=handles.img;
axes(handles.axes2);
h=fspecial('prewitt');
g2=imfilter(handles.img,h);
g3=imadd(g2,handles.img);
end
imshow(g3);
handles.img=g3;
guidata(hObject,handles);

 %锐化laplacian滤波
function radiobutton8_Callback(hObject, eventdata, handles)
global T
if handles.filebig==0
msgbox('请输入函数图像！','error');
return;
else
 axes(handles.axes2);
 imshow(handles.img);
T=handles.img;
axes(handles.axes2);
h=fspecial('laplacian');
g2=imfilter(handles.img,h);
g3=imadd(g2,handles.img);
imshow(g3);
end
handlse.img=g3;
guidata(hObject,handles);
