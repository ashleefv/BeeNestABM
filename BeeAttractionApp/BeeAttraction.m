function varargout = BeeAttraction(varargin)
% BEEATTRACTION MATLAB code for BeeAttraction.fig
%      BEEATTRACTION, by itself, creates a new BEEATTRACTION or raises the existing
%      singleton*.
%
%      H = BEEATTRACTION returns the handle to a new BEEATTRACTION or the handle to
%      the existing singleton*.
%
%      BEEATTRACTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BEEATTRACTION.M with the given input arguments.
%
%      BEEATTRACTION('Property','Value',...) creates a new BEEATTRACTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BeeAttraction_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BeeAttraction_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BeeAttraction

% Last Modified by GUIDE v2.5 14-Nov-2017 23:41:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BeeAttraction_OpeningFcn, ...
                   'gui_OutputFcn',  @BeeAttraction_OutputFcn, ...
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


% --- Executes just before BeeAttraction is made visible.
function BeeAttraction_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BeeAttraction (see VARARGIN)
% default parameters
handles.attraction = 0.0;
handles.totalTimePoints = 150;
handles.output = hObject;
addpath(genpath('beefiles'))
axes(handles.axes1);
bee = imread('beefiles/Beeimage_app.png');
imshow(bee);
% Update handles structuree
guidata(hObject,handles)

% UIWAIT makes BeeAttraction wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BeeAttraction_OutputFcn(hObject, eventdata, handles) 
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
% calls the simulation folder.
axes(handles.axes3);
simulation_attraction_app(hObject, eventdata, handles)




function perin_Callback(hObject, eventdata, handles)
% hObject    handle to perin (see GCBO) %percent in
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
in = str2double(get(handles.perin,'String'));
if in <= 100
    in = in/100;
    handles.attraction = in;
    %Saves the handles
    guidata(hObject,handles)
else
    errordlg('Please pick an actual percentage. Simulation invalid.','Input Error')
%Saves the handles
guidata(hObject,handles)

% Hints: get(hObject,'String') returns contents of perin as text
%        str2double(get(hObject,'String')) returns contents of perin as a double
end

% --- Executes during object creation, after setting all properties.
function perin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to perin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Determine the selected data set.
str = get(hObject,'String');
val = get(hObject,'Value');
%Sets the length of time

switch str{val};
    case 'Short' 
        handles.totalTimePoints = 150
    case 'Medium'
        handles.totalTimePoints = 300
    case 'Long'
        handles.totalTimePoints = 1500
end

%Saves the handles
guidata(hObject,handles)
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
