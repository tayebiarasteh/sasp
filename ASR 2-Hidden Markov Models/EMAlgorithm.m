function varargout = EMAlgorithm(varargin)
% EMALGORITHM MATLAB code for EMAlgorithm.fig
%      EMALGORITHM, by itself, creates a new EMALGORITHM or raises the existing
%      singleton*.
%
%      H = EMALGORITHM returns the handle to a new EMALGORITHM or the handle to
%      the existing singleton*.
%
%      EMALGORITHM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EMALGORITHM.M with the given input arguments.
%
%      EMALGORITHM('Property','Value',...) creates a new EMALGORITHM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EMAlgorithm_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EMAlgorithm_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EMAlgorithm

% Last Modified by GUIDE v2.5 03-Sep-2015 15:25:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @EMAlgorithm_OpeningFcn, ...
    'gui_OutputFcn',  @EMAlgorithm_OutputFcn, ...
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


% --- Executes just before EMAlgorithm is made visible.
function EMAlgorithm_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EMAlgorithm (see VARARGIN)

% Choose default command line output for EMAlgorithm
handles.output = hObject;
hObject.Position(4) = 1100;

%%%%%%%% Set global parameters here:

%delay between steps of EM-algorithm in seconds
% handles.delay = 0.5;
handles.delay = 0.01;

%resolution of support vector for elipses
handles.resolution = 0.01;

%rounding for history axis
handles.scale_factor_history = 10;

%testsamples used for calculation of RR
handles.runns = 1000;

% %testsamples used for calculation of ROC
% handles.runns_roc = 1000;
% 
% %points used to calculate ROC
% handles.datapoint_roc = 10;

%%%%%%%% end of global parameters

%add HMM functions and EM functions
addpath('HMM');
addpath('EM_files');

%init mu of both GMM words, will be updated in each EM step
handles.mu1=[1,1;2,2;-1,4;3.5,2];
handles.mu2=handles.mu1;

%init sigma of both GMM words, will be updated in each EM step
sigma1=[1,0;0,1];
sigma2=[1,0;0,1];
sigma3=[1,0;0,1];
sigma4=[1,0;0,1];

handles.sigma1=sigma1;
handles.sigma1(:,:,2)=sigma2;
handles.sigma1(:,:,3)=sigma3;
handles.sigma1(:,:,4)=sigma4;

handles.sigma2=sigma1;
handles.sigma2(:,:,2)=sigma2;
handles.sigma2(:,:,3)=sigma3;
handles.sigma2(:,:,4)=sigma4;

%copy initial mu+sigma to GUI
handles.s1init.Data=handles.sigma1(:,:,1);
handles.s2init.Data=handles.sigma1(:,:,2);
handles.s3init.Data=handles.sigma1(:,:,3);
handles.s4init.Data=handles.sigma1(:,:,4);

handles.m1init.Data=handles.mu1(1,:);
handles.m2init.Data=handles.mu1(2,:);
handles.m3init.Data=handles.mu1(3,:);
handles.m4init.Data=handles.mu1(4,:);

%create space to save all past mu and initialze
handles.mu_history1=[];
handles.mu_history2=[];
handles.mu_history1(:,:,1)=handles.mu1;
handles.mu_history2(:,:,1)=handles.mu2;

%create space to save all past sigma an initialize
handles.sigma_history1=[];
handles.sigma_history2=[];
handles.sigma_history1(:,:,:,1)=handles.sigma1;
handles.sigma_history2(:,:,:,1)=handles.sigma2;

%create space to save word feature vectors
handles.data1 = [];
handles.data2 = [];

%create GMM for each word, mixing proportin is default (equal)
handles.pdf1=gmdistribution(handles.mu1,handles.sigma1);
handles.pdf2=gmdistribution(handles.mu2,handles.sigma2);

%Variable to save if GMM is trained
handles.trained1 = 0;
handles.trained2 = 0;

%Space to save testword for classifier
handles.testword=[];
handles.testwordnum = 0;

%init axes
config_data_axes(handles);
config_EM_history_axes(handles,20);
config_mixing_axes(handles);
config_classifier_axes(handles);

%init configuration variables
set(handles.popup_num_gauss,'Value',3);

%write to console:
'initialized'

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EMAlgorithm wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = EMAlgorithm_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in button_sim_w1.
function button_sim_w1_Callback(hObject, eventdata, handles)
%reset GMM
handles = reset_GMM(handles,1);

%simulate featurevector for word 1
handles.data1=simulate_word1();

%plot featurevectors
axes(handles.TS_ax1);
cla;
plot(handles.data1(1,:),handles.data1(2,:),'b.','color','c');

%scale axes
config_data_axes(handles)

% Update handles structure
guidata(hObject, handles);
% hObject    handle to button_sim_w1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in button_train_w1.
function button_train_w1_Callback(hObject, eventdata, handles)

%test if parameters are correct
if(~check_input_EM(handles))
    return
end
if(isempty(handles.data1))
    errordlg(['training data is not initialized']);
    return
end

%disable input to window while programm is running
block_window(true,handles);

%number of EM-steps to do
stepnum = str2num(get(handles.step_num,'String'));

%write on console
['doing ',num2str(stepnum),' EM-steps']

handles = reset_GMM(handles,1);

%plot initial values
plot_em_distribution(handles,1);
plot_change(handles,1,stepnum);
pause(handles.delay);


for i=1:stepnum
    %show current step in editable field for stepnumber
    set(handles.step_num,'String',num2str(i));
    
    %perform one EM-step
    handles = step_em(handles,1);
    
    %plot distibution + change
    plot_em_distribution(handles,1);
    plot_change(handles,1,stepnum);
    pause(handles.delay);
end

%enable input to window
block_window(false,handles);

handles.trained1 = 1;

%
% % Update handles structure
guidata(hObject, handles);



% hObject    handle to button_train_w1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in button_train_w2.
function button_train_w2_Callback(hObject, eventdata, handles)

%test if parameters are correct
if(~check_input_EM(handles))
    return
end
if(isempty(handles.data2))
    errordlg(['training data is not initialized']);
    return
end

%disable input to window while programm is running
block_window(true,handles);

%number of EM-steps to do
stepnum = str2num(get(handles.step_num,'String'));

%write on console
['doing ',num2str(stepnum),' EM-steps']

handles = reset_GMM(handles,2);

%plot initial values
plot_em_distribution(handles,2);
plot_change(handles,2,stepnum);
pause(handles.delay);


for i=1:stepnum
    %show current step in editable field for stepnumber
    set(handles.step_num,'String',num2str(i));
    
    %perform one EM-step
    handles = step_em(handles,2);
    
    %plot distibution + change
    plot_em_distribution(handles,2);
    plot_change(handles,2,stepnum);
    pause(handles.delay);
end

%enable input to window
block_window(false,handles);

handles.trained2 = 1;

%
% % Update handles structure
guidata(hObject, handles);

% hObject    handle to button_train_w2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in button_sim_w2.
function button_sim_w2_Callback(hObject, eventdata, handles)
%reset GMM
handles = reset_GMM(handles,2);

%simulate featurevector for word 2
handles.data2=simulate_word2();

%plot featurevectors
axes(handles.TS_ax2);
cla;
plot(handles.data2(1,:),handles.data2(2,:),'b.','color','c');

%scale axes
config_data_axes(handles)

% Update handles structure
guidata(hObject, handles);
% hObject    handle to button_sim_w2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function step_num_Callback(hObject, eventdata, handles)
% hObject    handle to step_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Catch errors when entering the stepnumber (only numbers between 1 and 100
%are valid)
stepnum = str2num(get(handles.step_num,'String'));
if(isempty(stepnum))
    set(handles.step_num,'String','20');
    errordlg('only numbers between 1 and 50 are valid');
elseif(stepnum > 100)
    set(handles.step_num,'String','100');
    errordlg('only up to 100 steps are possible');
elseif(stepnum <1)
    set(handles.step_num,'String','1');
    errordlg('only possitive stepnumbers are possible');
end

% Hints: get(hObject,'String') returns contents of step_num as text
%        str2double(get(hObject,'String')) returns contents of step_num as a double


% --- Executes during object creation, after setting all properties.
function step_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to step_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes when entered data in editable cell(s) in s1init.
function s1init_CellEditCallback(hObject, eventdata, handles)

%set both covariances to get a symetric matrix
data = get(hObject,'Data');
data(eventdata.Indices(2),eventdata.Indices(1)) = str2num(eventdata.EditData);
data = set(hObject,'Data',data);

check_input_EM(handles);

% hObject    handle to s1init (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when entered data in editable cell(s) in s2init.
function s2init_CellEditCallback(hObject, eventdata, handles)

%set both covariances to get a symetric matrix
data = get(hObject,'Data');
data(eventdata.Indices(2),eventdata.Indices(1)) = str2num(eventdata.EditData);
data = set(hObject,'Data',data);

check_input_EM(handles);
% hObject    handle to s2init (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when entered data in editable cell(s) in s3init.
function s3init_CellEditCallback(hObject, eventdata, handles)

%set both covariances to get a symetric matrix
data = get(hObject,'Data');
data(eventdata.Indices(2),eventdata.Indices(1)) = str2num(eventdata.EditData);
data = set(hObject,'Data',data);

check_input_EM(handles);
% hObject    handle to s3init (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when entered data in editable cell(s) in s4init.
function s4init_CellEditCallback(hObject, eventdata, handles)

%set both covariances to get a symetric matrix
data = get(hObject,'Data');
data(eventdata.Indices(2),eventdata.Indices(1)) = str2num(eventdata.EditData);
data = set(hObject,'Data',data);

check_input_EM(handles);
% hObject    handle to s4init (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when entered data in editable cell(s) in m1init.
function m1init_CellEditCallback(hObject, eventdata, handles)

check_input_EM(handles);
% hObject    handle to m1init (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

% --- Executes when entered data in editable cell(s) in m2init.
function m2init_CellEditCallback(hObject, eventdata, handles)

check_input_EM(handles);
% hObject    handle to m2init (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when entered data in editable cell(s) in m3init.
function m3init_CellEditCallback(hObject, eventdata, handles)

check_input_EM(handles);
% hObject    handle to m3init (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when entered data in editable cell(s) in m4init.
function m4init_CellEditCallback(hObject, eventdata, handles)

check_input_EM(handles);
% hObject    handle to m4init (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in test_word1.
function test_word1_Callback(hObject, eventdata, handles)

%save type of testword and create testword
handles.testwordnum = 1;
handles.testword = simulate_test_sample(1);

%plot featurevectors:
plot_classifier_axes(handles,1);
plot_classifier_axes(handles,2);

%write word to GUI
set(handles.testword_num,'string','word1');

%plot classifier bounds if aviable
plot_classifier();

%
% % Update handles structure
guidata(hObject, handles);

% hObject    handle to test_word1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in test_word2.
function test_word2_Callback(hObject, eventdata, handles)

%save type of testword and create testword
handles.testwordnum = 2;
handles.testword = simulate_test_sample(2);

%plot featurevectors:
plot_classifier_axes(handles,1);
plot_classifier_axes(handles,2);

%write word to GUI
set(handles.testword_num,'string','word2');

%plot classifier bounds if aviable
plot_classifier();

%
% % Update handles structure
guidata(hObject, handles);
% hObject    handle to test_word2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in classify.
function classify_Callback(hObject, eventdata, handles)


%test if a testword has been created
if(handles.testwordnum ==0)
    errordlg('No testword has been created');
    return
end

%check in priors are valid
if(check_input_classifier(handles)==0)
    return
end

%read priors
prior1 = str2num(get(handles.prior_w1,'string'));
prior2 = str2num(get(handles.prior_w2,'string'));


%start classification
[class,logLikeliehood] = bayes_classifier(handles.pdf1,handles.pdf2,[prior1,prior2],handles.testword);

set(handles.result_classification,'string',['word',num2str(class)]);
set(handles.log_like1,'string',num2str(logLikeliehood(1)));
set(handles.log_like2,'string',num2str(logLikeliehood(2)));


% hObject    handle to classify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function prior_w1_Callback(hObject, eventdata, handles)

if(str2num(get(handles.prior_w1,'string'))<0||str2num(get(handles.prior_w2,'string'))<0||str2num(get(handles.prior_w1,'string'))>1||str2num(get(handles.prior_w2,'string'))>1)
    set(handles.prior_w1,'string',num2str(1-str2num(get(handles.prior_w2,'string'))));
    errordlg('The priors must be numbers between 0 and 1');
    return
end

set(handles.prior_w2,'string',num2str(1-str2num(get(handles.prior_w1,'string'))));

% hObject    handle to prior_w1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prior_w1 as text
%        str2double(get(hObject,'String')) returns contents of prior_w1 as a double



function prior_w2_Callback(hObject, eventdata, handles)

if(str2num(get(handles.prior_w1,'string'))<0||str2num(get(handles.prior_w2,'string'))<0||str2num(get(handles.prior_w1,'string'))>1||str2num(get(handles.prior_w2,'string'))>1)
    set(handles.prior_w2,'string',num2str(1-str2num(get(handles.prior_w1,'string'))));
    errordlg('The priors must be numbers between 0 and 1');
    return
end

set(handles.prior_w1,'string',num2str(1-str2num(get(handles.prior_w2,'string'))));
% hObject    handle to prior_w2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prior_w2 as text
%        str2double(get(hObject,'String')) returns contents of prior_w2 as a double




%%%%%%%%%%%%%%functions not directly connectet to GUI

function handles = step_em(handles,wordnum)
%This Function performes a single Step of the EM-clustering algorithm

if(wordnum==1)
    %call function in EM_files
    handles.pdf1=EM_step(handles.data1,handles.pdf1);
    
    %save current sigma and mu
    handles.sigma1=handles.pdf1.Sigma;
    handles.mu1=handles.pdf1.mu;
    
    %add sigma and mu to saved history
    handles.mu_history1(:,:,end+1)=handles.mu1;
    handles.sigma_history1(:,:,:,end+1)=handles.sigma1;
else
    %call function in EM_files
    handles.pdf2=EM_step(handles.data2,handles.pdf2);
    
    %save current sigma and mu
    handles.sigma2=handles.pdf2.Sigma;
    handles.mu2=handles.pdf2.mu;
    
    %add sigma and mu to saved history
    handles.mu_history2(:,:,end+1)=handles.mu2;
    handles.sigma_history2(:,:,:,end+1)=handles.sigma2;
end

function plot_em_distribution(handles,wordnum)
%plots the training featurevectors, the current Gaussian distribution as an
%ellipsis with the principal components as semiaxes and all past
%meanvectors of the Gaussdistributions

%select and clear axes
axes(handles.(['TS_ax',num2str(wordnum)]));
ax = gca;
cla(ax);
hold on

%plot training featurevectors
plot(handles.(['data',num2str(wordnum)])(1,:),handles.(['data',num2str(wordnum)])(2,:),'b.','color','c')

%support vector to plot elipses
t=-pi:handles.resolution:pi;

%perform plot for each Gaussdistribution
for j=1:handles.(['pdf',num2str(wordnum)]).NumComponents
    
    %read sigma and mu
    mu=handles.(['mu',num2str(wordnum)])(j,:);
    sigma=handles.(['sigma',num2str(wordnum)])(:,(2*j-1):(2*j));
    
    %to plot the ellipsis we need the sizes of the semiaxes, which are the
    %principal components of the Gaussdistribution. The eigenvectors
    %describe the direction of the axes, the eigenvalues the length.
    %perform PCA
    [eigvect,eigval]=eig(sigma);
    [tmp,maxpos]=max(mu);
    
    %calculate rotation angle of elipsis
    phi=atan(eigvect(2,1)/eigvect(1,1));
    
    %calculate unrotated ellipsis with zero mean
    x=sqrt(eigval(1,1))*cos(t);
    x(2,:)=sqrt(eigval(2,2))*sin(t);
    
    %calculate transformation matrix to rotate elipses perform
    %transformation
    trans=[cos(phi),-sin(phi);sin(phi),cos(phi)];
    x=trans*x;
    
    %shift elipsis to the mean of the distribution
    x(1,:)=x(1,:)+mu(1);
    x(2,:)=x(2,:)+mu(2);
    
    %plot the elipsis
    plot(x(1,:),x(2,:),'LineWidth',2)
    
    %reset color to previous one to plot the past meanvectors in the same
    %color
    ax.ColorOrderIndex=ax.ColorOrderIndex-1;
    
    %plot all past meanvectors of the distribution
    plot(reshape(handles.(['mu_history',num2str(wordnum)])(j,1,:),1,[]),reshape(handles.(['mu_history',num2str(wordnum)])(j,2,:),1,[]),'Marker','x','LineWidth',1)
end

%scale and label axes
config_data_axes(handles);
hold off

axes(handles.(['Mixing_ax',num2str(wordnum)]));
ax = gca;
cla(ax);

mixing = handles.(['pdf',num2str(wordnum)]).ComponentProportion;
mixing(2,:)=0;

h = barh([0,10],mixing,'stacked');

cor = get(gca,'colororder');

for i=1:size(mixing,2)
    set(h(i), 'FaceColor', cor(i,:));
end

config_mixing_axes(handles);

function plot_change(handles,wordnum,stepnum)
%this function plots the change of the mean and covariance matrix,
%stepnumber is needed to scale the axes

%select axes to plot Mu change
axes(handles.Mu_ax);
ax = gca;
%clear axes
cla(ax);
hold on;

%read change in mu for all past spets
mu_change = handles.(['mu_history',num2str(wordnum)])(:,:,2:end) - handles.(['mu_history',num2str(wordnum)])(:,:,1:end-1);

%calculate euclidean norm of mean change for each timestep
mu_change = mu_change.^2;
mu_change_mean = [];
for i=1:size(mu_change,3)
    mu_change_mean = [mu_change_mean , mu_change(:,:,i)*[1;1]];
end
mu_change_mean = mu_change_mean.^0.5;

%plot change
plot(1:size(mu_change_mean,2),mu_change_mean','LineWidth',2)



%select axes to plot Sigma change
axes(handles.Sigma_ax);
ax = gca;
%clear axes
cla(ax);
hold on;

%read change in Sigma for all past spets
sigma_change = handles.(['sigma_history',num2str(wordnum)])(:,:,:,2:end) - handles.(['sigma_history',num2str(wordnum)])(:,:,:,1:end-1);

%calculate Frobenius norm of Sigma change for each timestep
sigma_change = sigma_change.^2;
sigma_change_mean = [];
for i=1:size(sigma_change,4)
    tmp = [];
    for j=1:size(sigma_change,3)
        tmp(j,1) = [1,1]*sigma_change(:,:,j,i)*[1;1];
    end
    sigma_change_mean = [sigma_change_mean , tmp];
end

mu_change_mean = mu_change_mean.^0.5;

%plot change
plot(1:size(sigma_change_mean,2),sigma_change_mean','LineWidth',2)

%set up axes display
config_EM_history_axes(handles,stepnum,max(max(mu_change_mean(:,2:end))),max(max(sigma_change_mean(:,2:end))));

hold off

function plot_classifier_axes(handles,axnum)

axes(handles.(['BC_ax',num2str(axnum)]));
ax = gca;
cla(ax);
hold on

plot(handles.testword(1,:),handles.testword(2,:),'b.','color','b');

if(handles.(['trained',num2str(axnum)])==0)
    return
end

%support vector to plot elipses
t=-pi:handles.resolution:pi;

%perform plot for each Gaussdistribution
for j=1:handles.(['pdf',num2str(axnum)]).NumComponents
    
    %read sigma and mu
    mu=handles.(['mu',num2str(axnum)])(j,:);
    sigma=handles.(['sigma',num2str(axnum)])(:,(2*j-1):(2*j));
    
    %to plot the ellipsis we need the sizes of the semiaxes, which are the
    %principal components of the Gaussdistribution. The eigenvectors
    %describe the direction of the axes, the eigenvalues the length.
    %perform PCA
    [eigvect,eigval]=eig(sigma);
    [tmp,maxpos]=max(mu);
    
    %calculate rotation angle of elipsis
    phi=atan(eigvect(2,1)/eigvect(1,1));
    
    %calculate unrotated ellipsis with zero mean
    x=sqrt(eigval(1,1))*cos(t);
    x(2,:)=sqrt(eigval(2,2))*sin(t);
    
    %calculate transformation matrix to rotate elipses perform
    %transformation
    trans=[cos(phi),-sin(phi);sin(phi),cos(phi)];
    x=trans*x;
    
    %shift elipsis to the mean of the distribution
    x(1,:)=x(1,:)+mu(1);
    x(2,:)=x(2,:)+mu(2);
    
    %plot the elipsis
    plot(x(1,:),x(2,:),'LineWidth',2)
end

%scale and label axes
config_data_axes(handles);
hold off

function config_data_axes(handles)

%init axes
axes(handles.TS_ax1);
ylim([0 5.5]);
xlim([-3 4]);
xlabel('OFFSET');
ylabel('POWER');
grid on;

axes(handles.TS_ax2);
ylim([0 5.5]);
xlim([-3 4]);
xlabel('OFFSET');
ylabel('POWER');
grid on;

function handles = reset_GMM(handles,wordnum)
%this function resets the GMM for the selected word to its initial values

handles.(['trained',num2str(wordnum)]) = 0;


num_gauss = get(handles.popup_num_gauss,'value')

%read initial values from GUI
for i=1:num_gauss;
    handles.(['sigma',num2str(wordnum)])(:,:,i)=handles.(['s',num2str(i),'init']).Data;
    handles.(['mu',num2str(wordnum)])(i,:)=handles.(['m',num2str(i),'init']).Data;
end

%create new GMM, mixing proportion is default(equal)
handles.(['pdf',num2str(wordnum)])=gmdistribution(handles.(['mu',num2str(wordnum)])(1:num_gauss,:),handles.(['sigma',num2str(wordnum)])(:,:,1:num_gauss));

%save current Sigma and Mu
handles.(['sigma',num2str(wordnum)])=handles.(['pdf',num2str(wordnum)]).Sigma;
handles.(['mu',num2str(wordnum)])=handles.(['pdf',num2str(wordnum)]).mu;

%Save current Sigma and Mu to plot all past values
handles.(['mu_history',num2str(wordnum)])=[];
handles.(['mu_history',num2str(wordnum)])(:,:,1)=handles.(['mu',num2str(wordnum)]);

handles.(['sigma_history',num2str(wordnum)])=[];
handles.(['sigma_history',num2str(wordnum)])(:,:,:,1)=handles.(['sigma',num2str(wordnum)]);

function block_window(setting,handles)
%turns turn off or on the GUI-input

if(setting)
    %turn off the interface
    InterfaceObj=findobj(handles.figure1,'Enable','on');
    set(InterfaceObj,'Enable','off');
else
    %turn on the interface
    InterfaceObj=findobj(handles.figure1,'Enable','off');
    set(InterfaceObj,'Enable','on');
end

function config_EM_history_axes(handles,stepnum, varargin)
%scale and lable axes to plot past Mu and Sigma

scale = handles.scale_factor_history;

axes(handles.Sigma_ax);
xlim([1 stepnum]);
xlabel('number of iteration');
ylabel('Frobenius norm of change in \Sigma');
grid on;
if(nargin == 4)
    if(~isempty(varargin{2}))
        ylim([0,ceil(varargin{2}*scale)/scale]);
    end
end

axes(handles.Mu_ax);
xlim([1 stepnum]);
xlabel('number of iteration');
ylabel('Euclidean norm of change in \mu');
grid on;
if(nargin == 4)
    if(~isempty(varargin{1}))
        ylim([0,ceil(varargin{1}*scale)/scale]);
    end
end

function result = check_input_EM(handles)
%this function checks if all needed user input is correct
result = 1;
for i = 1:4
    %read initial Sigma and Mu
    sigma = handles.(['s',num2str(i),'init']).Data;
    mu = handles.(['m',num2str(i),'init']).Data;
    
    %test if Sigma and mu only contain numbers
    if(~isnumeric(sigma)||~isnumeric(mu)||max(max(isnan(sigma)))||max(isnan(mu)))
        result = 0;
        errordlg(['Initial values for distribution ',num2str(i),' contain nonnumeric elements']);
    end
    
    %test if variance is positive
    if(sigma(1,1)<0||sigma(2,2)<0)
        result = 0;
        errordlg(['Initial values for distribution ',num2str(i),' contain negative variance values']);
    end
    
    %test if covariance matrix is symetric
    if(sigma(1,2)~=sigma(2,1))
        result = 0;
        errordlg(['Covariance matrix of distribution ',num2str(i),' must be symetric']);
    end
    
    %test if covariance matrix is positive definit
    if((sigma(1,1)*sigma(2,2)-sigma(1,2)^2)<=0)
        result = 0;
        errordlg({['Covariance matrix of distribution ',num2str(i),' must be postive semidefinit:'],['The absolute value of the covariance must be smaller than: ',num2str(sigma(1,1)*sigma(2,2))]});
    end
end

function config_mixing_axes(handles)
axes(handles.Mixing_ax1);
xlim([0,1]);
ylim([0,1]);
xlabel('mixing weights');
grid off;
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);

axes(handles.Mixing_ax2);
xlim([0,1]);
ylim([0,1]);
xlabel('mixing weights');
grid off;
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);

function plot_classifier()

function result = check_input_classifier(handles)
result = 1;
if(str2num(get(handles.prior_w1,'string'))<0||str2num(get(handles.prior_w2,'string'))<0||str2num(get(handles.prior_w1,'string'))>1||str2num(get(handles.prior_w2,'string'))>1)
    result = 0;
    errordlg('The priors must be numbers between 0 and 1');
end

if(abs(str2num(get(handles.prior_w1,'string'))+str2num(get(handles.prior_w2,'string'))-1)>0.0001)
    result = 0;
    errordlg('The summ of the priors must be 1');
end

%test if both GMM are trained
if(handles.trained1==0||handles.trained2==0)
    result = 0;
    errordlg('The GMM must be trained first');
end

function config_classifier_axes(handles)

%init axes
axes(handles.BC_ax1);
ylim([0 5.5]);
xlim([-3 4]);
xlabel('OFFSET');
ylabel('POWER');
grid on;

%init axes
axes(handles.BC_ax2);
ylim([0 5.5]);
xlim([-3 4]);
xlabel('OFFSET');
ylabel('POWER');
grid on;

% %init axes
% axes(handles.ROC_ax);
% ylim([0 1]);
% xlim([0 1]);
% xlabel('recall word1');
% ylabel('recall word2');
% set(gca,'XDir','reverse');
% grid on;

function recall = measure_recall(handles,wordnum)
if(check_input_classifier(handles)==0)
    recall = -1;
return
end

%read priors
prior1 = str2num(get(handles.prior_w1,'string'));
prior2 = str2num(get(handles.prior_w2,'string'));

%store errors
error_count = 0;

%waitbar
progbar = waitbar(0.0,['measuring recall for word',num2str(wordnum)]);
set(progbar,'WindowStyle','modal');

for i = 1:handles.runns/2

%start classification
[class,logLikeliehood] = bayes_classifier(handles.pdf1,handles.pdf2,[prior1,prior2],simulate_test_sample(wordnum));
if(class ~= wordnum)
    error_count = error_count +1;
end
%update waitbar only if 1% changed
if(mod(i,round(handles.runns/200))==0)
progbar = waitbar(2*i/handles.runns);
end
end
%close waitbar
close;
%calculate recognition rate
recall = (handles.runns/2-error_count)/(handles.runns/2);

% --- Executes on button press in measure_rr.
function measure_rr_Callback(hObject, eventdata, handles)
recall1 = measure_recall(handles,1);
recall2 = measure_recall(handles,2);

%check for error
if(recall1==-1||recall2==-1)
    return
end

RR = str2num(get(handles.prior_w1,'string'))*recall1 + str2num(get(handles.prior_w2,'string'))*recall2;

set(handles.recognition_rate,'string',[num2str(RR*100),'%']);

%%%%%%%%%removed functions

% hObject    handle to measure_rr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in measure_roc.
% function measure_roc_Callback(hObject, eventdata, handles)
% 
% save_runns = handles.runns;
% handles.runns = handles.runns_roc;
% 
% % Update handles structure
% guidata(hObject, handles);
% 
% 
% recall = zeros(2,10);
% 
% 
% for i=0:handles.datapoint_roc
%     
% param = sqrt(abs(cos(pi*i/handles.datapoint_roc)))*sign(cos(pi*i/handles.datapoint_roc));
%     
% set(handles.prior_w1,'string',num2str(0.5 - 0.5*param));
% set(handles.prior_w2,'string',num2str(0.5 + 0.5*param));
%     
%     %%%%intuitive method bad scaling
% % set(handles.prior_w1,'string',num2str(i/handles.datapoint_roc));
% % set(handles.prior_w2,'string',num2str(1-i/handles.datapoint_roc));
%     
% recall(1,i+1) = measure_recall(handles,1);
% recall(2,i+1) = measure_recall(handles,2);
% 
% %check for error
% if(recall(1,i+1)==-1||recall(2,i+1)==-1)
%     return
% end
% end
% 
% handles.runns = save_runns;
% 
% recall
% 
% axes(handles.ROC_ax);
% cla(gca);
% hold on;
% plot(recall(1,:),recall(2,:))
% plot([0,1],[1,0]);
% config_classifier_axes(handles)
% 
% % Update handles structure
% guidata(hObject, handles);

% hObject    handle to measure_roc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% function plot_decision_region(handles,priors)
% 
% %classify test data
% [tmp1,tmp2,log]=bayes_classifier(handles.pdf1,handles.pdf2,priors,handles.testword);
% 
% %axis limits + calculation precission
% x=[-2,4];
% y=[0,4];
% z=[-2,1];
% precission = 0.05;
% 
% %create figure
% figure('units','normalized','outerposition',[0.1,0.1,0.9,0.9]);
% hold on
% 
% %calculate decision surface
% dec_surface = calculate_decision(x,y,precission,handles.pdf1,handles.pdf2,priors);
% 
% %plot decision surface
% surf(x(1):precission:x(2),y(1):precission:y(2),dec_surface,'LineWidth',0.01,'EdgeAlpha',0.4)
% %factor to use full colormap
% colormap_factor = max(max(abs(dec_surface)));
% caxis([-colormap_factor, colormap_factor]);
% %select colormap
% colormap(jet);
% 
% %Plot decision boundary
% contour3(x(1):precission:x(2),y(1):precission:y(2),dec_surface,[0,0],'LineWidth',2,'Color','black')
% 
% plot3(handles.testword(1,:),handles.testword(2,:),log,'b.','color','r');
% 
% %axis setup
% colorbar();
% xlim(x);
% ylim(y);
% zlim(z);
% view(145,45);
% axis equal;

% % --- Executes on button press in show_decision.
% function show_decision_Callback(hObject, eventdata, handles)
% %test if a testword has been created
% if(handles.testwordnum ==0)
%     errordlg('No testword has been created');
%     return
% end
% 
% %check in priors are valid
% if(check_input_classifier(handles)==0)
%     return
% end
% 
% %read priors
% prior1 = str2num(get(handles.prior_w1,'string'));
% prior2 = str2num(get(handles.prior_w2,'string'));
% 
% plot_decision_region(handles,[prior1,prior2]);
% % hObject    handle to show_decision (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
