function varargout = IsolatedWordRec(varargin)
% ISOLATEDWORDREC M-file for IsolatedWordRec.fig
%      ISOLATEDWORDREC, by itself, creates hmm1_a new ISOLATEDWORDREC or raises the existing
%      singleton*.
%
%      H = ISOLATEDWORDREC returns the handle to hmm1_a new ISOLATEDWORDREC or the handle to
%      the existing singleton*.
%
%      ISOLATEDWORDREC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ISOLATEDWORDREC.M with the given input arguments.
%
%      ISOLATEDWORDREC('Property','Value',...) creates hmm1_a new ISOLATEDWORDREC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before IsolatedWordRec_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to IsolatedWordRec_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help IsolatedWordRec

% Last Modified by GUIDE v2.5 23-Sep-2011 10:14:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IsolatedWordRec_OpeningFcn, ...
                   'gui_OutputFcn',  @IsolatedWordRec_OutputFcn, ...
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


% --- Executes just before IsolatedWordRec is made visible.
function IsolatedWordRec_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in hmm1_a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to IsolatedWordRec (see VARARGIN)

% Choose default command line output for IsolatedWordRec
handles.output = hObject;


set(handles.disp_A1,'Visible','Off');
set(handles.disp_pi1,'Visible','Off');
set(handles.disp_A1,'Position',[51 115 273 66]);
set(handles.disp_pi1,'Position',[51 40 273 24]);
handles.Active_plotHMM1 = 0;
handles.Active_trainTS1 = 0;

set(handles.disp_A2,'Visible','Off');
set(handles.disp_pi2,'Visible','Off');
set(handles.disp_A2,'Position',[51 115 273 66]);
set(handles.disp_pi2,'Position',[51 40 273 24]);
handles.Active_plotHMM2 = 0;
handles.Active_trainTS2 = 0;

set(handles.bg_LL1,'Visible','Off')
set(handles.bg_LL2,'Visible','Off')
set(handles.HMM1pdf1,'Visible','Off')
set(handles.HMM1pdf2,'Visible','Off')
set(handles.HMM1pdf3,'Visible','Off')
set(handles.HMM2pdf1,'Visible','Off')
set(handles.HMM2pdf2,'Visible','Off')
set(handles.HMM2pdf3,'Visible','Off')

handles.available_HMM1 = 0;
handles.available_HMM2 = 0;
handles.available_TestSeq = 0;

% Update handles structure
guidata(hObject, handles);

path('./HMM/',path)

% UIWAIT makes IsolatedWordRec wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = IsolatedWordRec_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in hmm1_a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in simTS1.
function simTS1_Callback(hObject, eventdata, handles)
% hObject    handle to simTS1 (see GCBO)
% eventdata  reserved - to be defined in hmm1_a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles = guihandles(hObject); 
handles.Active_plotHMM1 = 0;
handles.available_HMM1 = 0;
set(handles.disp_A1,'Visible','Off');
set(handles.disp_pi1,'Visible','Off');
set(handles.bg_LL1,'Visible','Off')
set(handles.bg_LL2,'Visible','Off')
axes(handles.HMM1pdf1); cla;
axes(handles.HMM1pdf2); cla;
axes(handles.HMM1pdf3); cla;

% State transition and state priors
A=[.8 .1 .1;
   .1 .8 .1;
   .1 .1 .8];

Pi=[ 1 0 0];

% create 'nrecord' records, each with 'nsteps' steps.
nrecord=10;
nsteps=400; % time steps per record
N=16;       % length of each segment
NFEAT=2;

[x,istart,nsamp]=hmm_maketestdata1(Pi,A,nrecord,nsteps,N,NFEAT);
axes(handles.axesTS1)
plot(handles.axesTS1,x(1,:),x(2,:),'b.');
set(handles.axesTS1,'xlim',[-3 4])
set(handles.axesTS1,'ylim',[0 4.5])
xlabel('OFFSET');
ylabel('POWER');
handles.xref1 = x;
handles.istart = istart;
handles.nsamp = nsamp;

handles.Active_trainTS1 = 1;

guidata(hObject, handles);



% --- Executes on button press in trainTS1.
function trainTS1_Callback(hObject, eventdata, handles)
% hObject    handle to trainTS1 (see GCBO)
% eventdata  reserved - to be defined in hmm1_a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.Active_trainTS1
	notFinished = 1;
	while notFinished
		try
			set(handles.bg_LL1,'Visible','Off')
			set(handles.bg_LL2,'Visible','Off')
			x = handles.xref1;
			istart = handles.istart;
			nsamp = handles.nsamp;
			% Initialize HMM parameters
			names={'OFFSET','POWER'};
			min_std=[.1 .1];
			NSTATES=3;
			NMODE=10;
			parm=init_hmm(x,NSTATES,NMODE,names,min_std);
			
			% Train
			NIT=100; % number of iterations
			[log_pdf_val, parm] = hmm2_reest(parm, x, istart, nsamp, NIT);
			
			set(handles.disp_A1,'Data',round(parm.A*1000)/1000);
			set(handles.disp_A1,'Visible','On');
			set(handles.disp_pi1,'Data',(round(parm.Pi*1000)/1000)');
			set(handles.disp_pi1,'Visible','On');
			
			handles.Active_plotHMM1 = 1;
			handles.parm1 = parm;
			handles.available_HMM1 = 1;
			
			figure(111);close(111)
			p = parm.pdf(1);axes(handles.HMM1pdf1);gmix_view2(p,x,1,2,0,200);axis off;
			p = parm.pdf(2);axes(handles.HMM1pdf2);gmix_view2(p,x,1,2,0,200);axis off;
			p = parm.pdf(3);axes(handles.HMM1pdf3);gmix_view2(p,x,1,2,0,200);axis off;
			
			guidata(hObject, handles);
			notFinished = 0;
		catch err
			notFinished = 1;
		end
	end
else
	msgbox('No training data available. Please simulate a training set first.','Missing data','warn','modal')
end

% --- Executes on button press in simTS2.
function simTS2_Callback(hObject, eventdata, handles)
% hObject    handle to simTS2 (see GCBO)
% eventdata  reserved - to be defined in hmm1_a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Active_plotHMM2 = 0;
handles.available_HMM2 = 0;
set(handles.disp_A2,'Visible','Off');
set(handles.disp_pi2,'Visible','Off');
set(handles.bg_LL1,'Visible','Off')
set(handles.bg_LL2,'Visible','Off')
axes(handles.HMM2pdf1); cla;
axes(handles.HMM2pdf2); cla;
axes(handles.HMM2pdf3); cla;

% State transition and state priors
A=[.7 .2 .1;
   .2 .7 .1;
   .1 .1 .8];

Pi=[ 1 0 0];

% create 'nrecord' records, each with 'nsteps' steps.
nrecord=10;
nsteps=400; % time steps per record
N=16;       % length of each segment
NFEAT=2;

[x,istart,nsamp]=hmm_maketestdata2(Pi,A,nrecord,nsteps,N,NFEAT);
axes(handles.axesTS2)
plot(handles.axesTS2,x(1,:),x(2,:),'b.');
set(handles.axesTS2,'xlim',[-3 4])
set(handles.axesTS2,'ylim',[0 4.5])
xlabel('OFFSET');
ylabel('POWER');
handles.xref2 = x;
handles.istart = istart;
handles.nsamp = nsamp;

handles.Active_trainTS2 = 1;

guidata(hObject, handles);

% --- Executes on button press in showHMM1.
function showHMM1_Callback(hObject, eventdata, handles)
% hObject    handle to showHMM1 (see GCBO)
% eventdata  reserved - to be defined in hmm1_a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotHMM1(handles);

function plotHMM1(handles)
if(handles.Active_plotHMM1)
% 	dialog('position',[100 100 1200 400]);
	figure('position',[100 100 1200 800]);
	parm = handles.parm1;
	x = handles.xref1;
	hmm2_view(parm,x,1,2);
	subplot(2,3,2)
	plot(x(1,:),x(2,:),'b.');
	xlabel('OFFSET');
	ylabel('POWER');
	title('Training set');
else
	msgbox('No HMM data available. Train an HMM before displaying it.','Missing data','warn','modal')
end


% --- Executes on button press in trainTS2.
function trainTS2_Callback(hObject, eventdata, handles)
% hObject    handle to trainTS2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.Active_trainTS2
	notFinished = 1;
	while notFinished
		try
			set(handles.bg_LL1,'Visible','Off')
			set(handles.bg_LL2,'Visible','Off')
			x = handles.xref2;
			istart = handles.istart;
			nsamp = handles.nsamp;
			% Initialize HMM parameters
			names={'OFFSET','POWER'};
			min_std=[.1 .1];
			NSTATES=3;
			NMODE=10;
			parm=init_hmm(x,NSTATES,NMODE,names,min_std);
			
			% Train
			NIT=100; % number of iterations
			[log_pdf_val, parm] = hmm2_reest(parm, x, istart, nsamp, NIT);
			
			set(handles.disp_A2,'Data',round(parm.A*1000)/1000);
			set(handles.disp_A2,'Visible','On');
			set(handles.disp_pi2,'Data',(round(parm.Pi*1000)/1000)');
			set(handles.disp_pi2,'Visible','On');
			
			handles.Active_plotHMM2 = 1;
			handles.parm2 = parm;
			handles.available_HMM2 = 1;
			
			
			p = parm.pdf(1);axes(handles.HMM2pdf1);gmix_view2(p,x,1,2,0,200);axis off;
			p = parm.pdf(2);axes(handles.HMM2pdf2);gmix_view2(p,x,1,2,0,200);axis off;
			p = parm.pdf(3);axes(handles.HMM2pdf3);gmix_view2(p,x,1,2,0,200);axis off;
			
			guidata(hObject, handles);
			notFinished = 0;
		catch err
			notFinished = 1;
		end
	end
else
	msgbox('No training data available. Please simulate a training set first.','Missing data','warn','modal')
end

% --- Executes on button press in showHMM2.
function showHMM2_Callback(hObject, eventdata, handles)
% hObject    handle to showHMM2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotHMM2(handles);

function plotHMM2(handles)
if(handles.Active_plotHMM2)
% 	dialog('position',[100 100 1200 400]);
	figure('position',[100 100 1200 800]);
	parm = handles.parm2;
	x = handles.xref2;
	hmm2_view(parm,x,1,2);
	subplot(2,3,2)
	plot(x(1,:),x(2,:),'b.');
	xlabel('OFFSET');
	ylabel('POWER');
	title('Training set');
else
	msgbox('No HMM data available. Train an HMM before displaying it.','Missing data','warn','modal')
end


% --- Executes on button press in createTestSeq1.
function createTestSeq1_Callback(hObject, eventdata, handles)
% hObject    handle to createTestSeq1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.bg_LL1,'Visible','Off')
set(handles.bg_LL2,'Visible','Off')
set(handles.disp_LL1,'String','')
set(handles.disp_LL2,'String','')

testword = 1;
hmm_example_generate_testseq;
axes(handles.axesTestSeq)
plot(handles.axesTestSeq,x2(1,:),x2(2,:),'b*');
set(handles.axesTestSeq,'xlim',[-3 4])
set(handles.axesTestSeq,'ylim',[0 4.5])
xlabel('OFFSET');
ylabel('POWER');
title('Test sequence');
set(handles.disp_testWord,'String',['Word ' num2str(testword)]);
handles.x2 = x2;
handles.available_TestSeq = 1;
guidata(hObject, handles);


% --- Executes on button press in classify.
function classify_Callback(hObject, eventdata, handles)
% hObject    handle to classify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~handles.available_HMM1
	msgbox('HMM for word 1 not trained.','Missing data','warn','modal')
elseif ~handles.available_HMM2
	msgbox('HMM for word 2 not trained.','Missing data','warn','modal')
elseif ~handles.available_TestSeq
	msgbox('No test sequence available.','Missing data','warn','modal')
else
	x2 = handles.x2;
	parm = handles.parm1;
	hmm_example_classify_testseq;
	q1 = q;
	parm = handles.parm2;
	hmm_example_classify_testseq;
	q2 = q;
	set(handles.disp_LL1,'String',num2str(q1))
	set(handles.disp_LL2,'String',num2str(q2))
	if q1>q2
		set(handles.bg_LL1,'Visible','On')
	else
		set(handles.bg_LL2,'Visible','On')
	end
end


% --- Executes on button press in createTestSeq2.
function createTestSeq2_Callback(hObject, eventdata, handles)
% hObject    handle to createTestSeq2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.bg_LL1,'Visible','Off')
set(handles.bg_LL2,'Visible','Off')
set(handles.disp_LL1,'String','')
set(handles.disp_LL2,'String','')

testword = 2;
hmm_example_generate_testseq;
axes(handles.axesTestSeq)
plot(handles.axesTestSeq,x2(1,:),x2(2,:),'b*');
set(handles.axesTestSeq,'xlim',[-3 4])
set(handles.axesTestSeq,'ylim',[0 4.5])
xlabel('OFFSET');
ylabel('POWER');
title('Test sequence');
set(handles.disp_testWord,'String',['Word ' num2str(testword)]);
handles.x2 = x2;
handles.available_TestSeq = 1;
guidata(hObject, handles);
