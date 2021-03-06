function varargout = JHJM(varargin)
% JHJM MATLAB code for JHJM.fig
%      JHJM, by itself, creates a new JHJM or raises the existing
%      singleton*.
%
%      H = JHJM returns the handle to a new JHJM or the handle to
%      the existing singleton*.
%
%      JHJM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in JHJM.M with the given input arguments.
%
%      JHJM('Property','Value',...) creates a new JHJM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before JHJM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to JHJM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help JHJM

% Last Modified by GUIDE v2.5 23-Mar-2022 13:26:44

% Begin initialization code - DO NOT EDIT
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
% End initialization code - DO NOT EDIT


% --- Executes just before JHJM is made visible.
function JHJM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to JHJM (see VARARGIN)

% Choose default command line output for JHJM
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes JHJM wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = JHJM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global s  %??????????????????????????????????????????????????????
[filename,pathname,filterindex]=...
uigetfile({'*.*';'*.bmp';'*.tif';'*.png';'*.jpg';'*.jpeg'},'select picture');  %??????????????????   
str=[pathname filename];  %????????????+?????????
s=str;
handles.filebig=filterindex;
if filterindex==0
%msgbox('?????????????????????','error');
return
else   
im=imread(str);   %????????????   
end 
axes(handles.axes1);  %???????????????axes
imshow(im);  %????????????
handles.img=im;
guidata(hObject,handles);


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5
 global T
if handles.filebig==0 %?????????????????????????????????
%msgbox('????????????????????????','error');
return;
else
 axes(handles.axes2);
 imshow(handles.img);
T=handles.img;
axes(handles.axes2);
prompt={'????????????????????????'};
defans={'3'};
p=inputdlg(prompt,'input',1,defans);
p1=str2num(p{1});
h1=fspecial('average',[p1 p1]);
I=imfilter(handles.img,h1);%??????????????????
end
imshow(I);
handles.img=I;
guidata(hObject,handles);


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton6
global T
if handles.filebig==0
%msgbox('????????????????????????','error');
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


% --- Executes on button press in radiobutton7.
function radiobutton7_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton7
global T
if handles.filebig==0
%msgbox('????????????????????????','error');
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


% --- Executes on button press in radiobutton8.
function radiobutton8_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton8
global T
if handles.filebig==0
msgbox('????????????????????????','error');
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
