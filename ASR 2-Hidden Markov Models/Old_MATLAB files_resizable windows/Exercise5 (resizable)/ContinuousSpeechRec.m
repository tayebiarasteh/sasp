function varargout = ContinuousSpeechRec(varargin)
% CONTINUOUSSPEECHREC M-file for ContinuousSpeechRec.fig
%      CONTINUOUSSPEECHREC, by itself, creates hmm1_a new CONTINUOUSSPEECHREC or raises the existing
%      singleton*.
%
%      H = CONTINUOUSSPEECHREC returns the handle to hmm1_a new CONTINUOUSSPEECHREC or the handle to
%      the existing singleton*.
%
%      CONTINUOUSSPEECHREC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONTINUOUSSPEECHREC.M with the given input arguments.
%
%      CONTINUOUSSPEECHREC('Property','Value',...) creates hmm1_a new CONTINUOUSSPEECHREC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ContinuousSpeechRec_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ContinuousSpeechRec_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ContinuousSpeechRec

% Last Modified by GUIDE v2.5 19-Oct-2011 09:06:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ContinuousSpeechRec_OpeningFcn, ...
                   'gui_OutputFcn',  @ContinuousSpeechRec_OutputFcn, ...
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


% --- Executes just before ContinuousSpeechRec is made visible.
function ContinuousSpeechRec_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in hmm1_a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ContinuousSpeechRec (see VARARGIN)

% Choose default command line output for ContinuousSpeechRec
handles.output = hObject;

path('./HMM/',path)

load parm
handles.parm = parm;
x=[-1.2 2; 0.5 3];
handles.xDummy = x;
set(handles.disp_A1,'Position',[0.1 0.3 0.6 0.2]);
set(handles.disp_pi1,'Position',[0.1 0.1 0.6 0.1]);
set(handles.disp_A1,'Data',[.8 .1 .1;.1 .8 .1; .1 .1 .8]);
set(handles.disp_pi1,'Data',[1 0 0]);

p = parm.pdf(1);axes(handles.HMM1pdf1);gmix_view2(p,x,1,2,0,200);axis off;
p = parm.pdf(2);axes(handles.HMM1pdf2);gmix_view2(p,x,1,2,0,200);axis off;
p = parm.pdf(3);axes(handles.HMM1pdf3);gmix_view2(p,x,1,2,0,200);axis off;

handles.available_HMM1 = 0;
handles.available_TestSeq = 0;

% Update handles structure
guidata(hObject, handles);

path('/HOMES/smeier/Lehre/SPSA/WS1112/supplements/MatlabSolution/Exercise6/HMM/',path)

% UIWAIT makes ContinuousSpeechRec wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ContinuousSpeechRec_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in hmm1_a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in showHMM1.
function showHMM1_Callback(hObject, eventdata, handles)
% hObject    handle to showHMM1 (see GCBO)
% eventdata  reserved - to be defined in hmm1_a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure(11)
set(11,'position',[100 100 1000 350]);
parm = handles.parm;
N=parm.N;
x = handles.xDummy;
for i=1:N,
	subplot(1,N,i);
	p = parm.pdf(i);
	gmix_view2(p,x,1,2,0,200);

	title(sprintf('State %d',i));
end;




% --- Executes on button press in createTestSeq.
function createTestSeq_Callback(hObject, eventdata, handles)
% hObject    handle to createTestSeq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.dispViterbi)
cla
testword = 1;
hmm_example_generate_testseq;
axes(handles.axesSeqFeatures)
plot(x2');
% legend('OFFSET','POWER');
xlabel('Time index')
axes(handles.axesSeqStates)
plot(states); ylim([.7 3.3]); 
set(gca,'XTickLabel',[]) 
set(gca,'YTick',[1 2 3]) 

handles.available_TestSeq = 1;
handles.x2 = x2;
handles.states = states;
guidata(hObject, handles);


% --- Executes on button press in startViterbi.
function startViterbi_Callback(hObject, eventdata, handles)
% hObject    handle to startViterbi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.available_TestSeq
	x2 = handles.x2;
	states = handles.states;
	parm = handles.parm;
	log_probs = hmm2_get_probs(parm,x2);
	est_states=hmm2_viterb(parm,log_probs,handles);
else
	msgbox('No feature sequence available for decoding.','Missing data','warn','modal')
end
	


% --- Executes on button press in fast.
function fast_Callback(hObject, eventdata, handles)
% hObject    handle to fast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fast
