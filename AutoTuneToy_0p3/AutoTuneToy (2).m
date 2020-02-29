function varargout = AutoTuneToy(varargin)
% AUTOTUNETOY M-file for AutoTuneToy.fig
%      AUTOTUNETOY, by itself, creates a new AUTOTUNETOY or raises the existing
%      singleton*.
%
%      H = AUTOTUNETOY returns the handle to a new AUTOTUNETOY or the handle to
%      the existing singleton*.
%
%      AUTOTUNETOY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AUTOTUNETOY.M with the given input arguments.
%
%      AUTOTUNETOY('Property','Value',...) creates a new AUTOTUNETOY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AutoTuneToy_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AutoTuneToy_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AutoTuneToy

% Last Modified by GUIDE v2.5 09-Jan-2010 16:20:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AutoTuneToy_OpeningFcn, ...
                   'gui_OutputFcn',  @AutoTuneToy_OutputFcn, ...
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


% --- Executes just before AutoTuneToy is made visible.
function AutoTuneToy_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AutoTuneToy (see VARARGIN)

% Choose default command line output for AutoTuneToy
handles.output = hObject;

handles.version = '0.1';

set(handles.axes1,'XTickLabel',{});
set(handles.axes1,'YTickLabel',{});


% Axes2 is just for mouse click input
set(handles.axes2,'XTick',[]);
set(handles.axes2,'YTick',[]);
set(handles.axes2,'XTickLabel',{});
set(handles.axes2,'YTickLabel',{});
set(handles.axes2,'ButtonDownFcn','AutoTuneToy(''mouseclick'',gcbo,[],guidata(gcbo))');


% Set defaults for the program
handles = set_defaults(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AutoTuneToy wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AutoTuneToy_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function handles = set_defaults(handles)

% Recording options defaults
handles.record_options.Fs = 11025;              % Sampling frequency (Hz)
handles.record_options.nbits = 16;              % Bits of precision
handles.record_options.precision = 'single';    % Save as single precision vector
handles.record_options.nchans = 1;              % Single audio channel (mono)
handles.record_options.initial_trim = 0.1;      % Always trim off the initial 0.1 sec

% Scale defaults
handles.scale_options.scale = 'Cmajor';         % Scale for pitch correction
handles.scale_options.indices = [];
handles.scale_options.freqs = [];
handles.scale_options.notes = {};
handles.scale_options.fund_index = [];

% Pitch detection defaults
handles.pitch_options.Fmin = 50;             	% Min frequency to search
handles.pitch_options.Fmax = 600;            	% Max frequency to search
handles.pitch_options.block_length = 0.04;   	% Length of each chuck to analyze for frequency
handles.pitch_options.step_per_block = 2;       % What fraction of the block width to step for each pitch detect
handles.pitch_options.threshold_amp = 0.25;      % Threshold of autocorr to be called a peak
handles.pitch_options.harmonic_deviation = 0.2; % Max deviation from multiple fundamental to be considered a harmonic (0.2 = +/- 20%)
handles.pitch_options.threshold_ratio = 0.6;    % Max ratio in heights of autocorr allowed between fundamental and harmonics
handles.pitch_options.f0_target_sample_time = 0.08; % Over what time interval to sample previous freq's to determine next target f0
handles.pitch_options.target_tol = 0.3;         % Max tolerance from f0 target (0.2 = +/-20%)

handles.pitch_options.slider1 = [0.02:0.01:0.1];
handles.pitch_options.slider2 = [1:8];
handles.pitch_options.slider3 = [0:0.05:1];
handles.pitch_options.slider4 = [0:0.05:1];
handles.pitch_options.slider5 = [0:0.05:1];
handles.pitch_options.slider6 = [0:0.02:0.4];
handles.pitch_options.slider7 = [0:0.05:1];
handles.pitch_options.fields = {
    'block_length'
    'step_per_block'
    'threshold_amp'
    'harmonic_deviation'
    'threshold_ratio'
    'f0_target_sample_time'
    'target_tol'
    };

handles.pitch_compress_options.scale_factor = 0.5;
handles.pitch_compress_options.decay_time = 0;

handles.pitch_compress_options.slider1 = [0:0.10:2];
handles.pitch_compress_options.slider2 = [0:0.5:10];
handles.pitch_compress_options.fields = {
    'scale_factor'
    'decay_time'
    };

handles.pitch_snap_options.snap_delay = 0;

handles.pitch_snap_options.slider1 = [0:0.1:1];
handles.pitch_snap_options.fields = {
    'snap_delay'
    };

handles.vibrato_options.amplitude = 0.2;
handles.vibrato_options.frequency = 6;
handles.vibrato_options.ramp = 0.4;

handles.vibrato_options.slider1 = [0:0.1:1];
handles.vibrato_options.slider2 = [1:1:10];
handles.vibrato_options.slider3 = [0:0.2:2];
handles.vibrato_options.fields = {
    'amplitude'
    'frequency'
    'ramp'
    };
    
% Initialize objects for record and playback
handles.record_obj = [];
handles.play_obj = [];

% Initialize structure for saving/manipulating sound data
handles.sound.A = [];
handles.sound.A_corrected = [];
handles.sound.Fs = [];
handles.sound.t = [];
handles.sound.f0 = [];
handles.sound.f0_save = [];
handles.sound.f0_corrected = [];
handles.sound.f0_corrected_save = [];
handles.sound.scale_factor = [];
handles.sound.tcalc = [];
handles.sound.selected_points = logical([]);
handles.sound.selected_points_save = logical([]);

% Initialize status reporting
handles.status.isrecording = false;
handles.status.isplaying = false;
handles.status.isreset = false;
handles.status.ispitch = false;
handles.status.iscompress = false;
handles.status.issnap = false;
handles.status.isvibrato = false;
handles.status.isundo = false;
handles.status.isview = false;
handles.status.ismousezoom = false;
handles.status.issave = false;
handles.status.Xlim = [];
handles.status.Ylim = [];
handles.status.filename = '';


[handles.scale_options.indices,...
    handles.scale_options.freqs,...
    handles.scale_options.notes,...
    handles.scale_options.fund_index] = ...
    get_scale(handles.scale_options.scale);

handles.status.isreset = true;
handles = update_GUI(handles);

set(handles.original_button,'value',1);
set(handles.plot_selection,'value',2);



function handles = update_plot(handles)

val = get(handles.plot_selection,'value');
axes(handles.axes1);

if val == 1
    if ~isempty(handles.sound.t) && ~isempty(handles.sound.A)
        hplot = [];
        legends = {};
        if ~isempty(handles.sound.A_corrected)
            hplot(1) = plot(handles.axes1,handles.sound.t,handles.sound.A,'color',[0.7 0.7 0.7]);
        else
            hplot(1) = plot(handles.axes1,handles.sound.t,handles.sound.A,'color','b');
        end
        legends{1} = 'Original';
        if ~isempty(handles.sound.A_corrected)
            hold on
            hplot(2) = plot(handles.axes1,handles.sound.t,handles.sound.A_corrected,'color','b');
            legends{2} = 'Modified';
            hold off
        end
        if isempty(handles.status.Xlim)
            set(handles.axes1,'Xlim',[min(handles.sound.t) max(handles.sound.t)]);
        else
            set(handles.axes1,'Xlim',handles.status.Xlim);
        end
        dA = max(handles.sound.A)-min(handles.sound.A);
        set(handles.axes1,'Ylim',[min(handles.sound.A)-0.1*dA max(handles.sound.A)+0.1*dA]);
        %legend(hplot,legends);
    else
        cla
    end
elseif val == 2
    if ~isempty(handles.sound.tcalc) && ~isempty(handles.sound.f0)
        hplot = [];
        f0 = handles.sound.f0;
        f0(f0 < 1) = NaN;
        if all(isnan(f0))
            cla
            hwarn = warndlg('No pitch detected with current settings!','WARNING');
            waitfor(hwarn);
        else
            hplot(1) = semilogy(handles.axes1,handles.sound.tcalc,f0,'.-','color',[0.7 0.7 0.7]);
            %legends{1} = 'Original';
            if ~isempty(handles.sound.f0_corrected)
                hold on
                f02 = handles.sound.f0_corrected;
                f02(f02 < 1) = NaN;
                hplot(2) = semilogy(handles.axes1,handles.sound.tcalc,f02,'.-','color','b');
                %legends{2} = 'Modified';
                
                if any(handles.sound.selected_points)
                    f03 = f02;
                    f03(~handles.sound.selected_points) = NaN;
                    semilogy(handles.axes1,...
                        handles.sound.tcalc,f03,'.-','color','r');
                end
                hold off
            end
            grid on
%             factors = [1.2,10];
%             tmp = min(f0)./factors;
%             Ylim(1) = max(tmp(tmp <= min(f02)));
%             tmp = max(f0).*factors;
%             Ylim(2) = min(tmp(tmp >= max(f02)));
            Ylim = [min(f0)/1.2 max(f0)*1.2];
            set(handles.axes1,'Ylim',Ylim);
            set(handles.axes1,'Ytick',handles.scale_options.freqs);
            set(handles.axes1,'YtickLabel',handles.scale_options.notes);
            if isempty(handles.status.Xlim)
                set(handles.axes1,'Xlim',[min(handles.sound.t) max(handles.sound.t)]);
            else
                set(handles.axes1,'Xlim',handles.status.Xlim);
            end
            %legend(hplot,legends);
        end
    else
        cla
    end
end

% Invisible axes for getting mouse information
set(handles.axes2,'Xlim',get(handles.axes1,'Xlim'));
set(handles.axes2,'Ylim',get(handles.axes1,'Ylim'));
set(handles.axes2,'Yscale',get(handles.axes1,'Yscale'));
set(handles.axes2,'ButtonDownFcn','AutoTuneToy(''mouseclick'',gcbo,[],guidata(gcbo))');
set(handles.axes2,'Color','none');
axes(handles.axes2);

clear f0 f02 f03

function handles = update_GUI(handles)

if handles.status.isrecording
    set(handles.play_button,'enable','off');
    set(handles.record_button,'enable','off');
    set(handles.stop_button,'enable','on');
elseif handles.status.isplaying
    set(handles.play_button,'enable','off');
    set(handles.record_button,'enable','off');
    set(handles.stop_button,'enable','off');
else
    if ~isempty(handles.sound.A)
        set(handles.play_button,'enable','on');
    else
        set(handles.play_button,'enable','off');
    end
    set(handles.record_button,'enable','on');
    set(handles.stop_button,'enable','off');
end
    

% Set all apply buttons to disabled until there is a frequency vector
if isempty(handles.sound.f0_corrected) || handles.status.isplaying
    state = 'off';
else
    state = 'on';
end
% set(handles.apply_compress_button,'enable',state);
% set(handles.snap_up_button,'enable',state);
% set(handles.snap_down_button,'enable',state);
% set(handles.apply_vibrato_button,'enable',state);
% set(handles.delete_button,'enable',state);
% set(handles.join_button,'enable',state);    
% set(handles.zoom_mouse_button,'enable',state); 
% set(handles.zoom_in_button,'enable',state); 
% set(handles.zoom_out_button,'enable',state); 
% set(handles.zoom_full_button,'enable',state); 
% set(handles.pan_left_button,'enable',state); 
% set(handles.pan_right_button,'enable',state); 
% set(handles.undo_button,'enable',state); 
% set(handles.undo_all_button,'enable',state); 
    
% Update pitch detection panel
if handles.status.isreset || handles.status.ispitch
    fields = handles.pitch_options.fields;
    
    set(handles.pitch_text1,'String',...
        sprintf('%0.0f ms',1000*handles.pitch_options.(fields{1})));
    set(handles.pitch_text2,'String',...
        sprintf('%0.0f',handles.pitch_options.(fields{2})));
    set(handles.pitch_text3,'String',...
        sprintf('%0.0f',100*handles.pitch_options.(fields{3})));
    set(handles.pitch_text4,'String',...
        sprintf('%0.0f %%',100*handles.pitch_options.(fields{4})));
    set(handles.pitch_text5,'String',...
        sprintf('%0.0f',100*handles.pitch_options.(fields{5})));
    set(handles.pitch_text6,'String',...
        sprintf('%0.0f ms',1000*handles.pitch_options.(fields{6})));
    set(handles.pitch_text7,'String',...
        sprintf('%0.0f %%',100*handles.pitch_options.(fields{7})));
    handles.status.ispitch = false;
end
if handles.status.isreset
    for N = 1:7
        all_vals = handles.pitch_options.(sprintf('slider%i',N));
        val = (handles.pitch_options.(fields{N}) - all_vals(1)) / ...
            (all_vals(end) - all_vals(1));
        set(handles.(sprintf('pitch_slider%i',N)),'value',val);
        set(handles.(sprintf('pitch_slider%i',N)),'sliderstep',...
            [1 1].*(1/(length(all_vals)-1)));
    end
end

% Update pitch compression panel
if handles.status.isreset || handles.status.iscompress
    fields = handles.pitch_compress_options.fields;
    
    set(handles.compress_text1,'String',...
        sprintf('%0.0f %%',100*handles.pitch_compress_options.(fields{1})));
    set(handles.compress_text2,'String',...
        sprintf('%0.1f s',handles.pitch_compress_options.(fields{2})));
    handles.status.iscompress = false;
end
if handles.status.isreset
    for N = 1:2
        all_vals = handles.pitch_compress_options.(sprintf('slider%i',N));
        val = (handles.pitch_compress_options.(fields{N}) - all_vals(1)) / ...
            (all_vals(end) - all_vals(1));
        set(handles.(sprintf('compress_slider%i',N)),'value',val);
        set(handles.(sprintf('compress_slider%i',N)),'sliderstep',...
            [1 1].*(1/(length(all_vals)-1)));
    end
end
 
% Update pitch snap panel
if handles.status.isreset || handles.status.issnap
    fields = handles.pitch_snap_options.fields;
    
    set(handles.snap_text1,'String',...
        sprintf('%0.1f s',handles.pitch_snap_options.(fields{1})));
    handles.status.issnap = false;
end
if handles.status.isreset
    for N = 1:1
        all_vals = handles.pitch_snap_options.(sprintf('slider%i',N));
        val = (handles.pitch_snap_options.(fields{N}) - all_vals(1)) / ...
            (all_vals(end) - all_vals(1));
        set(handles.(sprintf('snap_slider%i',N)),'value',val);
        set(handles.(sprintf('snap_slider%i',N)),'sliderstep',...
            [1 1].*(1/(length(all_vals)-1)));
    end
end

% Update vibrato panel
if handles.status.isreset || handles.status.isvibrato
    fields = handles.vibrato_options.fields;
    
    set(handles.vibrato_text1,'String',...
        sprintf('%0.1f',handles.vibrato_options.(fields{1})));
    set(handles.vibrato_text2,'String',...
        sprintf('%0.0f Hz',handles.vibrato_options.(fields{2})));
    set(handles.vibrato_text3,'String',...
        sprintf('%0.1f s',handles.vibrato_options.(fields{3})));
    handles.status.isvibrato = false;
end
if handles.status.isreset
    for N = 1:3
        all_vals = handles.vibrato_options.(sprintf('slider%i',N));
        val = (handles.vibrato_options.(fields{N}) - all_vals(1)) / ...
            (all_vals(end) - all_vals(1));
        set(handles.(sprintf('vibrato_slider%i',N)),'value',val);
        set(handles.(sprintf('vibrato_slider%i',N)),'sliderstep',...
            [1 1].*(1/(length(all_vals)-1)));
    end
end

% Fill in Scale Selection menu
if handles.status.isreset
    [indices,freqs,notes,fund_index] = get_scale();
    set(handles.scale_selection_list,'String',[indices,'CUSTOM']);
    I = find(strcmp(indices,handles.scale_options.scale));
    if ~isempty(I)
        set(handles.scale_selection_list,'Value',I);
    end
end

% Handle undo stuff
if isempty(handles.sound.f0_corrected_save);
    set(handles.undo_button,'enable','off');
    set(handles.undo_all_button,'enable','off');
else
    set(handles.undo_button,'enable','on');
    set(handles.undo_all_button,'enable','on');
end
handles.status.isundo = false;

% Handle view/zoom stuff
if handles.status.isview || handles.status.isreset
    if isempty(handles.sound.A)
        onoff = 'off';
    else
        onoff = 'on';
    end
    set(handles.plot_selection,'enable',onoff);
    set(handles.zoom_mouse_button,'enable',onoff);
    set(handles.zoom_in_button,'enable',onoff);
    set(handles.zoom_full_button,'enable',onoff);
    set(handles.pan_right_button,'enable',onoff);
    set(handles.pan_left_button,'enable',onoff);
    
    if isempty(handles.status.Xlim) && isempty(handles.status.Ylim)
        set(handles.zoom_out_button,'enable','off');
        set(handles.zoom_full_button,'enable','off');
    else
        set(handles.zoom_out_button,'enable','on');
        set(handles.zoom_full_button,'enable','on');
    end
    
    handles.status.isview = false;
end

if handles.status.isreset || handles.status.issave
    if isempty(handles.status.filename)
        set(handles.figure1,'Name','AutoTune Toy [untitled]');
    else
        [pathstr, name, ext, versn] = fileparts(handles.status.filename);
        set(handles.figure1,'Name',sprintf('AutoTune Toy [%s]',name))
    end
    handles.status.issave = false;
end

handles.status.isreset = false;    
 
function handles = update_pitch_sliders(handles,N)

fn = handles.pitch_options.fields{N};
slider_val = get(handles.(sprintf('pitch_slider%i',N)),'value');
all_vals = handles.pitch_options.(sprintf('slider%i',N));
new_val = slider_val*(all_vals(end)-all_vals(1))+all_vals(1);
[what,where] = min(abs(all_vals - new_val));
handles.pitch_options.(fn) = all_vals(where);

handles.status.ispitch = true;
handles = update_GUI(handles);

function handles = update_compress_sliders(handles,N)

fn = handles.pitch_compress_options.fields{N};
slider_val = get(handles.(sprintf('compress_slider%i',N)),'value');
all_vals = handles.pitch_compress_options.(sprintf('slider%i',N));
new_val = slider_val*(all_vals(end)-all_vals(1))+all_vals(1);
[what,where] = min(abs(all_vals - new_val));
handles.pitch_compress_options.(fn) = all_vals(where);

handles.status.iscompress = true;
handles = update_GUI(handles);

function handles = update_snap_sliders(handles,N)

fn = handles.pitch_snap_options.fields{N};
slider_val = get(handles.(sprintf('snap_slider%i',N)),'value');
all_vals = handles.pitch_snap_options.(sprintf('slider%i',N));
new_val = slider_val*(all_vals(end)-all_vals(1))+all_vals(1);
[what,where] = min(abs(all_vals - new_val));
handles.pitch_snap_options.(fn) = all_vals(where);

handles.status.issnap = true;
handles = update_GUI(handles);

function handles = update_vibrato_sliders(handles,N)

fn = handles.vibrato_options.fields{N};
slider_val = get(handles.(sprintf('vibrato_slider%i',N)),'value');
all_vals = handles.vibrato_options.(sprintf('slider%i',N));
new_val = slider_val*(all_vals(end)-all_vals(1))+all_vals(1);
[what,where] = min(abs(all_vals - new_val));
handles.vibrato_options.(fn) = all_vals(where);

handles.status.isvibrato = true;
handles = update_GUI(handles);    


% --- Executes on button press in 
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to record_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in record_button.
function record_button_Callback(hObject, eventdata, handles)
% hObject    handle to record_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    handles.record_obj = audiorecorder(handles.record_options.Fs,...
        handles.record_options.nbits,handles.record_options.nchans);
    record(handles.record_obj);

    handles.status.isrecording = true;

    handles = update_GUI(handles);
    guidata(hObject, handles);
catch
    h = errordlg('Error starting audio input device!  Check sound card and microphone!','ERROR');
    waitfor(h);
    return
end


% --- Executes on button press in stop_button.
function stop_button_Callback(hObject, eventdata, handles)
% hObject    handle to stop_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.status.isrecording
    stop(handles.record_obj);
    handles.status.isrecording = false;
    
    % Get the waveform
    handles.sound.A = getaudiodata(handles.record_obj,...
        handles.record_options.precision);
    
    % Save sampling frequency
    handles.sound.Fs = handles.record_options.Fs;
    
    % Trim initial samples
    Nstart = round(handles.record_options.initial_trim*...
        handles.record_options.Fs);
    handles.sound.A = handles.sound.A(Nstart:end);
    handles.sound.A_corrected = handles.sound.A;
    
    % Calculate time vector
    handles.sound.t = (0:length(handles.sound.A)-1)./handles.record_options.Fs;
    
    % Clear old frequency data
    handles.sound.f0 = [];
    handles.sound.f0_save = [];
    handles.sound.f0_corrected = [];
    handles.sound.f0_corrected_save = [];
    handles.sound.selected_points = logical([]);
    handles.sound.selected_points_save = logical([]);
    
    % Clear plot limits
    handles.status.Xlim = [];
    handles.status.Ylim = [];
    
    handles = run_pitch_detection(handles);
    handles.status.isview = true;
    handles = update_plot(handles);
    handles = update_GUI(handles);
    guidata(hObject, handles);
    
end



% --- Executes on button press in play_button.
function play_button_Callback(hObject, eventdata, handles)
% hObject    handle to play_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Decide which to play
val = get(handles.original_button,'val');
if val == 1
    if isempty(handles.sound.A)
        return
    end
    A = handles.sound.A;
else
    if isempty(handles.sound.A_corrected)
        handles = run_pitch_correction(handles);
    end
    A = handles.sound.A_corrected;
end



handles.status.isplaying = true;
handles = update_GUI(handles);
handles = update_plot(handles);

% Get the current plot limits
Ylim = get(handles.axes1,'Ylim');
if isempty(handles.status.Xlim)
    tmin = min(handles.sound.t);
    tmax = max(handles.sound.t);
else
    [what,Imin] = min(abs(handles.status.Xlim(1)-handles.sound.t));
    [what,Imax] = min(abs(handles.status.Xlim(2)-handles.sound.t));
    tmin = handles.sound.t(Imin);
    tmax = handles.sound.t(Imax);
    A = A(Imin:Imax);
end

% Draw the vertical line to scan across
axes(handles.axes1);
hline = line(tmin.*[1 1],[Ylim],'color','r');

pause(0.1);

% Start playing the wave file
% Note: Updated wavplay to audioplayer on 12/10/12 for cross-platform functionality
%wavplay(A,handles.sound.Fs,'async');
p = audioplayer(A, handles.sound.Fs); 
play(p);
tic;
while 1
    tnow = toc + tmin;
    if tnow > tmax
        break
    end
    % Scan the line
    set(hline,'XData',[tnow tnow]);
    pause(0.02);
end
set(hline,'XData',[tmin tmin]);

% Done
handles.status.isplaying = false;
handles = update_GUI(handles);

guidata(hObject, handles);

clear A;
axes(handles.axes2);


% --------------------------------------------------------------------
function record_options_menu_Callback(hObject, eventdata, handles)
% hObject    handle to record_options_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Fs = handles.record_options.Fs;
nbits = handles.record_options.nbits;
initial_trim = handles.record_options.initial_trim;

prompt = {
    'Sampling rate (1000 - 44100 Hz):'
    'Bits per sample (8, 16, or 24):'
    'Record start delay (0 - 1000 ms):'
    };
dlg_title = 'Record options';
num_lines = [1 35];
def = {
    num2str(Fs)
    num2str(nbits)
    num2str(initial_trim*1000)
    };
answer = inputdlg(prompt,dlg_title,num_lines,def);
if ~isempty(answer)
    if ~isempty(answer{1})
        Fs = round(str2num(answer{1}));
    end
    if ~isempty(answer{2})
        nbits = round(str2num(answer{2}));
    end
    if ~isempty(answer{3})
        initial_trim = round(str2num(answer{3}))/1000;
    end

    input_error = false;
    if Fs < 1000 || Fs > 44100
        Fs = 11025;
        input_error = true;
    end
    if ~ismember(nbits,[8,16,24]);
        nbits = 16;
        input_error = true;
    end
    if initial_trim < 0 || initial_trim > 1
        initial_trim = 0.1;
        input_error = true;
    end
    if input_error
        hwarn = warndlg('One or more values out of range, default values used!',...
            'WARNING');
        waitfor(hwarn);
    end
    
    if handles.record_options.Fs ~= Fs || ...
            handles.record_options.nbits ~= nbits || ...
            handles.record_options.initial_trim ~= initial_trim

%         if ~isempty(handles.sound.A)
%             button = questdlg('This will erase previous sound data.  Continue?',...
%                 'WARNING','Yes','No','Yes');
%             if strcmp(button,'No')
%                 return
%             end
%         end
        
        handles.record_options.Fs = Fs;
        handles.record_options.nbits = nbits;
        handles.record_options.initial_trim = initial_trim;

%         handles.sound.A = [];
%         handles.sound.A_corrected = [];
%         handles.sound.t = [];
%         handles.sound.f0 = [];
%         handles.sound.f0_corrected = [];
%         handles.sound.tcalc = [];
%         handles.sound.selected_points = logical([]);
        
        handles = update_plot(handles);

        guidata(hObject, handles);
    end
    
end

% --- Executes on selection change in plot_selection.
function plot_selection_Callback(hObject, eventdata, handles)
% hObject    handle to plot_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns plot_selection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plot_selection

if isempty(handles.sound.A_corrected)
    handles = run_pitch_correction(handles);
end
handles = update_plot(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function plot_selection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in original_button.
function original_button_Callback(hObject, eventdata, handles)
% hObject    handle to original_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of original_button

set(handles.original_button,'val',1);
set(handles.modified_button,'val',0);


% --- Executes on button press in modified_button.
function modified_button_Callback(hObject, eventdata, handles)
% hObject    handle to modified_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of modified_button

set(handles.original_button,'val',0);
set(handles.modified_button,'val',1);


% --- Executes on slider movement.
function pitch_slider1_Callback(hObject, eventdata, handles)
% hObject    handle to pitch_slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

N = 1;
handles = update_pitch_sliders(handles,N);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pitch_slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pitch_slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function pitch_slider2_Callback(hObject, eventdata, handles)
% hObject    handle to pitch_slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

N = 2;
handles = update_pitch_sliders(handles,N);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pitch_slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pitch_slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function pitch_slider3_Callback(hObject, eventdata, handles)
% hObject    handle to pitch_slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

N = 3;
handles = update_pitch_sliders(handles,N);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pitch_slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pitch_slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function pitch_slider4_Callback(hObject, eventdata, handles)
% hObject    handle to pitch_slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

N = 4;
handles = update_pitch_sliders(handles,N);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pitch_slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pitch_slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function pitch_slider5_Callback(hObject, eventdata, handles)
% hObject    handle to pitch_slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

N = 5;
handles = update_pitch_sliders(handles,N);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pitch_slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pitch_slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function pitch_slider6_Callback(hObject, eventdata, handles)
% hObject    handle to pitch_slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

N = 6;
handles = update_pitch_sliders(handles,N);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pitch_slider6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pitch_slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function pitch_slider7_Callback(hObject, eventdata, handles)
% hObject    handle to pitch_slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

N = 7;
handles = update_pitch_sliders(handles,N);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pitch_slider7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pitch_slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in apply_pitch_button.
function apply_pitch_button_Callback(hObject, eventdata, handles)
% hObject    handle to apply_pitch_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.sound.A)
    handles = run_pitch_detection(handles);
    handles = update_plot(handles);
    handles = update_GUI(handles);
    guidata(hObject, handles);
end


% --- Executes on slider movement.
function snap_slider1_Callback(hObject, eventdata, handles)
% hObject    handle to snap_slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

N = 1;
handles = update_snap_sliders(handles,N);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function snap_slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to snap_slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on button press in snap_down_button.
function snap_button_Callback(hObject, eventdata, handles)
% hObject    handle to snap_down_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.sound.f0_corrected) || ~any(handles.sound.selected_points) ||...
        isempty(handles.scale_options.freqs);
    return
end

snap_direction = get(gcbo,'String');

% Pull out the frequency and time of the selected points
I = find(handles.sound.selected_points);
sp = handles.sound.selected_points(I(1):I(end));
f0 = handles.sound.f0_corrected(I(1):I(end));
t = handles.sound.tcalc(I(1):I(end));
t = t-t(1);

% Come up with a envelope for compressing that goes from 1 down to the
% desired compression factor over the desired time interval
if length(t) < 2
    where = 1;
elseif handles.pitch_snap_options.snap_delay > t(2)
    [what,where] = min(abs(t - handles.pitch_snap_options.snap_delay));
else
    where = 1;
end

tmp = f0(where:end);
mean_f0 = 10^mean(log10(tmp(sp(where:end))));
if strcmp(snap_direction,'Down')
    Iu = find(handles.scale_options.freqs < 0.99*mean_f0);
    mean_new = max(handles.scale_options.freqs(Iu));
elseif strcmp(snap_direction,'Up')
    Iu = find(handles.scale_options.freqs > 1.01*mean_f0);
    mean_new = min(handles.scale_options.freqs(Iu));
end
factor = mean_new/mean_f0.*ones(size(t));
if where > 1
    factor(1:where) = linspace(1,factor(end),where);
end

% Scale the deviations around the mean in a log sense
f0 = f0.*factor;
handles.sound.f0_corrected_save(end+1,:) = handles.sound.f0_corrected;
handles.sound.f0_save(end+1,:) = handles.sound.f0;
handles.sound.selected_points_save(end+1,:) = handles.sound.selected_points;
handles.sound.f0_corrected(handles.sound.selected_points) = f0(sp);

handles.sound.A_corrected = [];

handles = update_plot(handles);
handles = update_GUI(handles);
guidata(hObject, handles);

clear t f0 what where factor

% --- Executes on slider movement.
function compress_slider1_Callback(hObject, eventdata, handles)
% hObject    handle to compress_slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

N = 1;
handles = update_compress_sliders(handles,N);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function compress_slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to compress_slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function compress_slider2_Callback(hObject, eventdata, handles)
% hObject    handle to compress_slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

N = 2;
handles = update_compress_sliders(handles,N);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function compress_slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to compress_slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in apply_compress_button.
function apply_compress_button_Callback(hObject, eventdata, handles)
% hObject    handle to apply_compress_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
                  
if isempty(handles.sound.f0_corrected) || ~any(handles.sound.selected_points)
    return
end

% Pull out the frequency and time of the selected points
f0 = handles.sound.f0_corrected(handles.sound.selected_points);
I = find(handles.sound.selected_points);
t = handles.sound.tcalc(I(1):I(end));
t = t-t(1);

% Come up with a envelope for compressing that goes from 1 down to the
% desired compression factor over the desired time interval
factor = ones(size(t)).*handles.pitch_compress_options.scale_factor;
if handles.pitch_compress_options.decay_time > t(2)
    [what,where] = min(abs(t - handles.pitch_compress_options.decay_time));
    factor(1:where) = linspace(1,handles.pitch_compress_options.scale_factor,where);
end
factor = factor(handles.sound.selected_points(I(1):I(end)));

% Scale the deviations around the mean in a log sense
f0 = 10.^(mean(log10(f0)) + factor.*(log10(f0)-mean(log10(f0))));
handles.sound.f0_corrected_save(end+1,:) = handles.sound.f0_corrected;
handles.sound.f0_save(end+1,:) = handles.sound.f0;
handles.sound.selected_points_save(end+1,:) = handles.sound.selected_points;
handles.sound.f0_corrected(handles.sound.selected_points) = f0;

handles.sound.A_corrected = [];

handles = update_plot(handles);
handles = update_GUI(handles);
guidata(hObject, handles);

clear t f0 what where factor



% --- Executes on slider movement.
function vibrato_slider1_Callback(hObject, eventdata, handles)
% hObject    handle to vibrato_slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

N = 1;
handles = update_vibrato_sliders(handles,N);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function vibrato_slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vibrato_slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function vibrato_slider2_Callback(hObject, eventdata, handles)
% hObject    handle to vibrato_slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

N = 2;
handles = update_vibrato_sliders(handles,N);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function vibrato_slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vibrato_slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function vibrato_slider3_Callback(hObject, eventdata, handles)
% hObject    handle to vibrato_slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

N = 3;
handles = update_vibrato_sliders(handles,N);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function vibrato_slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vibrato_slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in apply_vibrato_button.
function apply_vibrato_button_Callback(hObject, eventdata, handles)
% hObject    handle to apply_vibrato_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.sound.f0_corrected) || ~any(handles.sound.selected_points) ||...
        isempty(handles.scale_options.freqs);
    return
end

% Pull out the frequency and time of the selected points
I = find(handles.sound.selected_points);
sp = handles.sound.selected_points(I(1):I(end));
f0 = handles.sound.f0_corrected(I(1):I(end));
t = handles.sound.tcalc(I(1):I(end));
t = t-t(1);

% Come up with a envelope for compressing that goes from 1 down to the
% desired compression factor over the desired time interval
if length(t) < 2
    where = 1;
elseif handles.vibrato_options.ramp > t(2)
    [what,where] = min(abs(t - handles.vibrato_options.ramp));
else
    where = 1;
end

A = handles.vibrato_options.amplitude.*...
    sin(2*pi*handles.vibrato_options.frequency*t)/12;
if where > 1
    A(1:where) = A(1:where).*linspace(0,1,where);
end
factor = 2.^(A);

% Scale the deviations around the mean in a log sense
f0 = f0.*factor;
handles.sound.f0_corrected_save(end+1,:) = handles.sound.f0_corrected;
handles.sound.f0_save(end+1,:) = handles.sound.f0;
handles.sound.selected_points_save(end+1,:) = handles.sound.selected_points;
handles.sound.f0_corrected(handles.sound.selected_points) = f0(sp);

handles.sound.A_corrected = [];

handles = update_plot(handles);
handles = update_GUI(handles);
guidata(hObject, handles);

clear t f0 what where factor


% --- Executes on button press in delete_button.
function delete_button_Callback(hObject, eventdata, handles)
% hObject    handle to delete_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.sound.f0_corrected) || ~any(handles.sound.selected_points)
    return
end

% Save for undoing
handles.sound.f0_save(end+1,:) = handles.sound.f0;
handles.sound.f0_corrected_save(end+1,:) = handles.sound.f0_corrected;
handles.sound.selected_points_save(end+1,:) = handles.sound.selected_points;

% Pull out the frequency and time of the selected points
handles.sound.f0_corrected(handles.sound.selected_points) = 0;
handles.sound.f0(handles.sound.selected_points) = 0;
handles.sound.selected_points = false(size(handles.sound.f0));

handles.sound.A_corrected = [];

handles = update_plot(handles);
handles = update_GUI(handles);
guidata(hObject, handles);

clear t f0 what where factor

% --- Executes on button press in join_button.
function join_button_Callback(hObject, eventdata, handles)
% hObject    handle to join_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.sound.f0_corrected) || ~any(handles.sound.selected_points)
    return
end

% Save for undoing
handles.sound.f0_save(end+1,:) = handles.sound.f0;
handles.sound.f0_corrected_save(end+1,:) = handles.sound.f0_corrected;
handles.sound.selected_points_save(end+1,:) = handles.sound.selected_points;

% Pull out the frequency and time of the selected points
any_changes = false;
I = find(handles.sound.selected_points);
f0 = handles.sound.f0;
for n = 1:length(I)-1
    if I(n+1)-I(n) > 1
        if f0(I(n)) > 1 && f0(I(n+1)) > 1 && all(f0(I(n)+1:I(n+1)-1) < 1)
            f0(I(n):I(n+1)) = logspace(log10(f0(I(n))),log10(f0(I(n+1))),length(I(n):I(n+1)));
            handles.sound.selected_points(I(n):I(n+1)) = true;
            any_changes = true;
        end
    end
end
handles.sound.f0 = f0;

f0 = handles.sound.f0_corrected;
for n = 1:length(I)-1
    if I(n+1)-I(n) > 1
        if f0(I(n)) > 1 && f0(I(n+1)) > 1 && all(f0(I(n)+1:I(n+1)-1) < 1)
            f0(I(n):I(n+1)) = logspace(log10(f0(I(n))),log10(f0(I(n+1))),length(I(n):I(n+1)));
            any_changes = true;
        end
    end
end

if ~any_changes
    return
end

handles.sound.f0_corrected = f0;

handles.sound.A_corrected = [];

handles = update_plot(handles);
handles = update_GUI(handles);
guidata(hObject, handles);

clear f0

% --- Executes on button press in undo_button.
function undo_button_Callback(hObject, eventdata, handles)
% hObject    handle to undo_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.sound.f0_corrected_save)
    handles.sound.f0_corrected = handles.sound.f0_corrected_save(end,:);
    handles.sound.f0_corrected_save(end,:) = [];
    handles.sound.f0 = handles.sound.f0_save(end,:);
    handles.sound.f0_save(end,:) = [];
    handles.sound.selected_points = handles.sound.selected_points_save(end,:);
    handles.sound.selected_points_save(end,:) = [];
    
    handles.sound.A_corrected = [];
    
    handles = update_plot(handles);
    handles = update_GUI(handles);
    guidata(hObject, handles);
end


% --- Executes on button press in undo_all_button.
function undo_all_button_Callback(hObject, eventdata, handles)
% hObject    handle to undo_all_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.sound.f0_corrected_save)
    
    button = questdlg('Undo ALL changes?','Confirm undo','Yes','No','Yes');
    if strcmp(button,'No')
        return
    end
    
    handles.sound.f0_corrected = handles.sound.f0_corrected_save(1,:);
    handles.sound.f0_corrected_save = [];
    handles.sound.f0 = handles.sound.f0_save(1,:);
    handles.sound.f0_save = [];
    handles.sound.selected_points = handles.sound.selected_points_save(1,:);
    handles.sound.selected_points_save = logical([]);
    
    handles.sound.A_corrected = [];
    
    handles = update_plot(handles);
    handles = update_GUI(handles);
    guidata(hObject, handles);
end

% --- Executes on button press in redo_button.
function redo_button_Callback(hObject, eventdata, handles)
% hObject    handle to redo_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in zoom_in_button.
function zoom_in_button_Callback(hObject, eventdata, handles)
% hObject    handle to zoom_in_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.sound.A)
    return
end

if isempty(handles.status.Xlim)
    handles.status.Xlim = [min(handles.sound.t) max(handles.sound.t)];
end

dX = diff(handles.status.Xlim);
handles.status.Xlim(1) = handles.status.Xlim(1)+dX/10;
handles.status.Xlim(2) = handles.status.Xlim(2)-dX/10;
if diff(handles.status.Xlim) < diff(handles.sound.t(1:2))
    dt = diff(handles.sound.t(1:2));
    handles.status.Xlim = mean(handles.status.Xlim) + [-dt/2 dt/2];
end

handles.status.isview = true;

handles = update_plot(handles);
handles = update_GUI(handles);
guidata(hObject, handles);

% --- Executes on button press in zoom_mouse_button.
function zoom_mouse_button_Callback(hObject, eventdata, handles)
% hObject    handle to zoom_mouse_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.sound.A)
    return
end

handles.status.ismousezoom = true;
guidata(hObject, handles);

set(gcf,'Pointer','crosshair');


% --- Executes on button press in zoom_out_button.
function zoom_out_button_Callback(hObject, eventdata, handles)
% hObject    handle to zoom_out_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.sound.A)
    return
end

if isempty(handles.status.Xlim)
    handles.status.Xlim = [min(handles.sound.t) max(handles.sound.t)];
end

dX = diff(handles.status.Xlim);
handles.status.Xlim(1) = handles.status.Xlim(1)-dX/10;
handles.status.Xlim(2) = handles.status.Xlim(2)+dX/10;

% if (handles.status.Xlim(1) <= min(handles.sound.t)) && ...
%         (handles.status.Xlim(2) >= max(handles.sound.t))
%     handles.status.Xlim = [];
% else
%     if handles.status.Xlim(1) <= min(handles.sound.t)
%         handles.status.Xlim(1) = min(handles.sound.t);
%     end
%     if handles.status.Xlim(2) >= max(handles.sound.t)
%         handles.status.Xlim(2) = max(handles.sound.t);
%     end
% end

handles.status.isview = true;

handles = update_plot(handles);
handles = update_GUI(handles);
guidata(hObject, handles);

% --- Executes on button press in zoom_full_button.
function zoom_full_button_Callback(hObject, eventdata, handles)
% hObject    handle to zoom_full_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.sound.A)
    return
end

handles.status.Xlim = [];
handles.status.Ylim = [];
handles.status.isview = true;

handles = update_plot(handles);
handles = update_GUI(handles);
guidata(hObject, handles);

% --- Executes on button press in pan_left_button.
function pan_left_button_Callback(hObject, eventdata, handles)
% hObject    handle to pan_left_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.sound.A)
    return
end

if isempty(handles.status.Xlim)
    handles.status.Xlim = [min(handles.sound.t) max(handles.sound.t)];
end

dX = diff(handles.status.Xlim);
handles.status.Xlim(1) = handles.status.Xlim(1)-dX/10;
handles.status.Xlim(2) = handles.status.Xlim(2)-dX/10;

handles.status.isview = true;

handles = update_plot(handles);
handles = update_GUI(handles);
guidata(hObject, handles);

% --- Executes on button press in pan_right_button.
function pan_right_button_Callback(hObject, eventdata, handles)
% hObject    handle to pan_right_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.sound.A)
    return
end

if isempty(handles.status.Xlim)
    handles.status.Xlim = [min(handles.sound.t) max(handles.sound.t)];
end

dX = diff(handles.status.Xlim);
handles.status.Xlim(1) = handles.status.Xlim(1)+dX/10;
handles.status.Xlim(2) = handles.status.Xlim(2)+dX/10;

handles.status.isview = true;

handles = update_plot(handles);
handles = update_GUI(handles);
guidata(hObject, handles);

% --- Executes on selection change in scale_selection_list.
function scale_selection_list_Callback(hObject, eventdata, handles)
% hObject    handle to scale_selection_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns scale_selection_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from scale_selection_list


scale = get(handles.scale_selection_list,'String');
scale = scale{get(handles.scale_selection_list,'Value')};

if strcmp(scale,'CUSTOM')
    
else
    handles.scale_options.scale = scale;
    [handles.scale_options.indices,...
        handles.scale_options.freqs,...
        handles.scale_options.notes,...
        handles.scale_options.fund_index] = ...
        get_scale(handles.scale_options.scale);
end

handles = update_plot(handles);
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function scale_selection_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scale_selection_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data.sound = handles.sound;
data.status = handles.status;
data.record_options = handles.record_options;
data.scale_options = handles.scale_options;
data.pitch_options = handles.pitch_options;
data.pitch_compress_options = handles.pitch_compress_options;
data.pitch_snap_options = handles.pitch_snap_options;
data.vibrato_options = handles.vibrato_options;

if ~isempty(handles.status.filename)
    FilterSpec = handles.status.filename;
else
    FilterSpec = 'untitled.prj';
end
[FileName,PathName,FilterIndex] = uiputfile(FilterSpec,...
    'Select file to save');
if ischar(FileName)
    handles.status.filename = fullfile(PathName,FileName);
    data.status.filename = fullfile(PathName,FileName);
    save(fullfile(PathName,FileName), 'data','-mat');
end

handles.status.issave = true;
handles = update_GUI(handles);

guidata(hObject, handles);



% --------------------------------------------------------------------
function load_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

button = questdlg('Any changes since last save will be lost!  Continue?',...
    'Confirm load','Yes','No','Yes');
if strcmp(button,'No')
    return
end

if ~isempty(handles.status.filename)
    FilterSpec = handles.status.filename;
else
    FilterSpec = '*.prj';
end
[FileName,PathName,FilterIndex] = uigetfile(FilterSpec,...
    'Select file to load');
if ischar(FileName)
    load('-mat',fullfile(PathName,FileName));
    handles.status.filename = fullfile(PathName,FileName);


    handles.sound = data.sound;
    handles.status = data.status;
    handles.record_options = data.record_options;
    handles.scale_options = data.scale_options;
    handles.pitch_options = data.pitch_options;
    handles.pitch_compress_options = data.pitch_compress_options;
    handles.pitch_snap_options = data.pitch_snap_options;
    handles.vibrato_options = data.vibrato_options;



    handles.status.isreset = true;
    handles = update_GUI(handles);
    handles = update_plot(handles);

    guidata(hObject, handles);
end

% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function about_button_Callback(hObject, eventdata, handles)
% hObject    handle to about_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = msgbox({
    'AutoTune Toy'
    sprintf('Version %s',handles.version)
    'Carl Arft'
    'tfralrac@yahoo.com'
    },'About');
waitfor(h);


% --------------------------------------------------------------------
function help_button_Callback(hObject, eventdata, handles)
% hObject    handle to help_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


try
    % Changed "winopen" to "open" command on 12/10/12 (some install's could
    % not find README location)
    %winopen('README.pdf');
    open(fullfile(fileparts(which(mfilename)), 'README.pdf'));
catch
    h = errordlg({
        'Cannot open README.pdf.  Please locate and open'
        'this file manually for instructions on using the'
        'AutoTune Toy!!!'
        },'ERROR');
    waitfor(h)
end



function handles = run_pitch_detection(handles)
% hObject    handle to pitch_detect_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.sound.A)
    return
end

plot_on = false;

% Get variables out of the handles structure
block_length = handles.pitch_options.block_length;                  
step_per_block = handles.pitch_options.step_per_block;
Fmin = handles.pitch_options.Fmin;                            
Fmax = handles.pitch_options.Fmax;   
f0_target_sample_time = handles.pitch_options.f0_target_sample_time;
A = handles.sound.A;
t = handles.sound.t;
Fs = handles.sound.Fs;
target_tol = handles.pitch_options.target_tol;        
harmonic_deviation = handles.pitch_options.harmonic_deviation; 
threshold_ratio = handles.pitch_options.threshold_ratio;    
threshold_amp = handles.pitch_options.threshold_amp;     
 

block = step_per_block*round(block_length*Fs/step_per_block);     % Size of each block to find pitch
step = block/step_per_block;                         % Step (blocks are 4 times larger than "steps"

N = floor((length(A)-block)/step);      % Number of frequency computations

% Initialize variables for storing results
f0 = zeros(1,N);                        % Initialize vector for storing frequencies
tcalc = zeros(1,N);                     % The time at which that frequency calculation is valid
f0_target = [];                         % For keeping track of the target for the next calculation

f0_samples = floor...                   % Figure out how many f0 samples there will be
    (f0_target_sample_time/(step/Fs));

% Waitbar stuff
hwait = waitbar(0,'Determining pitch...');
set(hwait,'Name','Please Wait');
waitbar_count = 0;
waitbar_update = round(N/10);
pause(0.001);

%tic
I = 1;
for n = 1:N
    
    % Update waitbar
    if waitbar_count > waitbar_update
        waitbar(n/N,hwait);
        waitbar_count = 0;
    end
    waitbar_count = waitbar_count + 1;
    
    % Extract a block of the wav file
    Atemp = A(I:I+block-1);
    ttemp = t(I:I+block-1);
    I = I+step;
    
    % Do the autocorrelation to find the frequency
    [f0(n),acorr,Nshifts,Tshifts] = ...
        find_f0_timedomain2(Atemp,Fs,Fmin,Fmax,...
        threshold_amp,threshold_ratio,harmonic_deviation,...
        f0_target,target_tol);
    tcalc(n) = median(ttemp);
    
    if plot_on
        figure(2)
        subplot(2,1,1),
        plot(ttemp,Atemp)
        title(sprintf('Frequency: %0.1f Hz',f0(n)))
        subplot(2,1,2),
        plot(Tshifts,acorr)
        pause
    end
    
    % If the previous and current are invalid, then make the last one invalid
    % as well
    if n > 2
        if f0(n-2) < 1 && f0(n) < 1
            f0(n-1) = 0;
        end
    end
            
    % Look at the last few samples to determine the target frequency for
    % the next iteration
    if n >= f0_samples
        f0_chunk = f0(n-f0_samples+1:n);
        f0_target = f0_chunk(f0_chunk>0);
        if ~isempty(f0_target)
            f0_target = mean(f0_target);
        end
    end
    
    
end
%toc

if ishandle(hwait)
    close(hwait);
end

% Save calculation
handles.sound.A_corrected = [];
handles.sound.f0 = f0;
handles.sound.f0_save = [];
handles.sound.f0_corrected = f0;
handles.sound.f0_corrected_save= [];
handles.sound.tcalc = tcalc;
handles.sound.selected_points = false(size(f0));
handles.sound.selected_points_save= false(size(f0));
handles.sound.block = block;                  
handles.sound.step = step;
    
function handles = run_pitch_correction(handles)
% hObject    handle to pitch_detect_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.sound.A)
    return
end

plot_on = false;

A = handles.sound.A;
t = handles.sound.t;
f0 = handles.sound.f0;
f02 = handles.sound.f0_corrected;
Fs = handles.sound.Fs;
scale_factor = zeros(size(f0));
A_corrected = zeros(size(A));
block = handles.sound.block;                  
step = handles.sound.step;
N = length(f0);

max_acorr_shift = 0;
max_acorr_amp = 0;  

% Waitbar stuff
hwait = waitbar(0,'Correcting pitch...');
set(hwait,'Name','Please Wait');
waitbar_count = 0;
waitbar_update = round(N/10);
pause(0.001);

I = 1;
for n = 1:N
    
    % Update waitbar
    if waitbar_count > waitbar_update
        waitbar(n/N,hwait);
        waitbar_count = 0;
    end
    waitbar_count = waitbar_count + 1;

    % Pull out a chunk of the signal
    Atemp = A(I:I+block-1);
    ttemp = t(I:I+block-1);
       
    if f0(n) > 0 && ~isnan(f0(n))
        
        % Calculate the scale factor by which to shift the frequency
        scale_factor(n) = f02(n)/f0(n);
    
        % Interpolate
        tp = mean(ttemp) + (ttemp-mean(ttemp)).*scale_factor(n);
        Ainterp = interp1(ttemp,Atemp ,tp)';
        Ivalid = find(~isnan(Ainterp));
        Ainterp(isnan(Ainterp)) = 0;
        
        Nperiod = ceil(1/(f02(n)/Fs));
    else
        % No frequency shift
        scale_factor(n) = 1;
        Ainterp = Atemp;  
        Ivalid = find(~isnan(Ainterp));
        Nperiod = ceil(1/(handles.pitch_options.Fmin/Fs));
    end
    
    if n == 1
        A_corrected(I:I+block-1) = Ainterp;
        
    else
        
        % Pull out one period of the new waveform
        Achunk = Ainterp(Ivalid(1):Ivalid(1)+Nperiod-1);
        factor = sum(abs(Achunk))*2;
        
        if all(scale_factor(n-1:n) == 1)

        else

            % Start doing correlation
            max_acorr_amp = 0;
            max_acorr_shift = 0;

            for Nshift = -round(Nperiod/2)+ (1:Nperiod)

                % Calculate makeshift autocorrelation
                acorr = 1-sum(abs(Achunk - A_corrected((I:I+Nperiod-1) + Nshift + Ivalid(1))))./factor;

                if acorr > max_acorr_amp
                    max_acorr_amp = acorr;
                    max_acorr_shift = Nshift;
                end
            end
        end
        
        [what,where] = min(abs(Achunk - ...
            A_corrected((I:I+Nperiod-1) + max_acorr_shift + Ivalid(1))));
        
        if plot_on
            Asave = A_corrected((I:I+Nperiod-1) + max_acorr_shift + Ivalid(1));
        end
        
        A_corrected((I+where-1:I+length(Ivalid)-1) + max_acorr_shift + Ivalid(1)) = ...
            Ainterp(Ivalid(where:end));
        
        if plot_on
            tt = 1:length(Achunk);
            figure(1);
            plot(tt,Achunk,tt,Asave,tt,A_corrected((I:I+Nperiod-1) + max_acorr_shift + Ivalid(1)))
            title(sprintf('Acorr = %0.3f (Nshift = %i)',max_acorr_amp,max_acorr_shift))
            hold on
            plot(tt(where),Achunk(where),'*')
            hold off
            pause
        end

        

    end
    
    I = I+step;
    
end

handles.sound.A_corrected = A_corrected;
handles.sound.scale_factor = scale_factor;

clear f0 f02 scale factor A A_corrected

if ishandle(hwait)
    close(hwait);
end


% The following function makes nothing happen when you press a key
function nokeyresponse(hObject, eventdata, handles)

% The following function executes when the mouse is clicked on the plot


function mouseclick(hObject, eventdata, handles)


if isempty(handles.sound.f0)
    return
end

% keytype = get(gcf,'CurrentCharacter');
%
% if keytype == 's'

if handles.status.ismousezoom
    
    % Zoom with mouse 
    
    point1 = get(gca,'CurrentPoint');    % button down detected
    finalRect = rbbox;                  % return figure units
    point2 = get(gca,'CurrentPoint');    % button up detected
    point1 = point1(1,1:2);              % extract x and y
    point2 = point2(1,1:2);
    p1 = min(point1,point2);             % calculate locations
    offset = abs(point1-point2);         % and dimensions
    xminmax = [p1(1) p1(1)+offset(1)];
    yminmax = [p1(2) p1(2)+offset(2)];
    
    if diff(xminmax) < diff(handles.sound.t(1:2))
        dt = diff(handles.sound.t(1:2));
        xminmax = mean(xminmax) + [-dt/2 dt/2];
    end

    handles.status.Xlim = xminmax;
    handles.status.ismousezoom = false;
    handles.status.isview = true;
    
    handles = update_plot(handles);
    handles = update_GUI(handles);
    
    guidata(hObject, handles);
    set(gcf,'Pointer','arrow');
    
    return;

elseif get(handles.plot_selection,'value') ~= 1
    clicktype = get(gcf,'SelectionType');
    if strcmp(clicktype,'normal') || strcmp(clicktype,'alt')

        point1 = get(gca,'CurrentPoint');    % button down detected
        finalRect = rbbox;                  % return figure units
        point2 = get(gca,'CurrentPoint');    % button up detected
        point1 = point1(1,1:2);              % extract x and y
        point2 = point2(1,1:2);
        p1 = min(point1,point2);             % calculate locations
        offset = abs(point1-point2);         % and dimensions
        xminmax = [p1(1) p1(1)+offset(1)];
        yminmax = [p1(2) p1(2)+offset(2)];

        if strcmp(clicktype,'normal')
            handles.sound.selected_points(:) = false;
        end

        handles.sound.selected_points(...
            handles.sound.f0_corrected > yminmax(1) & ...
            handles.sound.f0_corrected < yminmax(2) & ...
            handles.sound.tcalc > xminmax(1) & ...
            handles.sound.tcalc < xminmax(2) & ...
            handles.sound.f0_corrected > 0) = true;

    elseif strcmp(clicktype,'extend')

        % Save for undoing
        handles.sound.f0_save(end+1,:) = handles.sound.f0;
        handles.sound.f0_corrected_save(end+1,:) = handles.sound.f0_corrected;
        handles.sound.selected_points_save(end+1,:) = handles.sound.selected_points;

        % Get initial click point
        point1 = get(gca,'CurrentPoint');    % button down detected
        point1 = point1(1,1:2);
        handles.point1 = point1;
        handles.single_point = false;

        % If no points have previously been selected, select the closest point
        if all(~handles.sound.selected_points)
            a = axis;
            tnorm = (handles.sound.tcalc - a(1))/(a(2)-a(1));
            fnorm = (handles.sound.f0_corrected - a(3))/(a(4)-a(3));
            xnorm = (point1(1) - a(1))/(a(2)-a(1));
            ynorm = (point1(2) - a(3))/(a(4)-a(3));
            [what,where] = min((tnorm-xnorm).^2+(fnorm-ynorm).^2);
            handles.sound.selected_points(where) = true;
            handles.single_point = true;
        end

        guidata(hObject, handles);
        pause(0.05);


        % Set function to track mouse position
        handles.tstart = now;
        guidata(hObject, handles);
        set(gcf,'WindowButtonMotionFcn','AutoTuneToy(''wbmcb'',gcbo,[],guidata(gcbo))');
        set(gcf,'WindowButtonUpFcn','AutoTuneToy(''wbucb'',gcbo,[],guidata(gcbo))');

        handles = update_GUI(handles);

    end
end

guidata(hObject, handles);
handles = update_plot(handles);


function wbmcb(hObject, eventdata, handles)
% Note: added in a timer so that the figure won't try to refresh too
% quickly and crash
if (now-handles.tstart)*24*60*60 > 0.08
    handles.tstart = now;
    
    cp = get(gca,'CurrentPoint');
    ydiff = cp(1,2)/handles.point1(2);
    handles.sound.f0_corrected(handles.sound.selected_points) = ...
        handles.sound.f0_corrected_save(end,handles.sound.selected_points).*ydiff;
    handles.sound.A_corrected = [];
    guidata(hObject, handles);
    handles = update_plot(handles);
else
    
end



function wbucb(hObject, eventdata, handles)
if handles.single_point
    handles.sound.selected_points(:) = false;

end
guidata(hObject, handles);
handles = update_plot(handles);
set(gcf,'WindowButtonMotionFcn','');
set(gcf,'WindowButtonUpFcn','');



function [f0,acorr,Nshifts,Tshifts] = find_f0_timedomain2(A,Fs,Fmin,Fmax,...
    threshold_amp,threshold_ratio,harmonic_deviation,target_f0,target_tol)

% This function calculates the fundamental frequency of a signal in the
% time domain

% OUTPUTS:
%
%   f0:         The interpolated fundamental freqency (Hz)
%   
%   acorr:      The vector of "autocorrelation" values
%
%   Nshifts:    The vector of shift indices associated with the values in
%               acorr
%
%   Tshifts:    The vector of time shifts associated with the values in
%               acorr
%
% INPUTS:
%   
%   A:          The input signal A(t), sampled at an interval Ts
%
%   Fs:         The sample frequency (Hz)
%
%   Fmin:       (Optional) The minimum expected frequency (Hz)
%
%   Fmax:       (Optional) The maximum expected freqency (Hz)
% 
%   threshold_amp: (Optional) The min autocorr amplitude considered peak
%   
%   threshold_ratio: (Optional) The min ratio between max autocorr peak and
%                       a given peak to be considered as a possible peak
%
%   harmonic deviation: (Optional) Tolerance for looking for another peak
%                       at half the frequency of the highest
%                       autocorrelation peak
%
%   target_f0:  (Optional) The probable target frequency
%
%   target_tol: (Optional) The max error between target_f0 and current freq
%                       (0.1 = 10%)
%

if nargin < 9 || isempty(target_tol)
    target_tol = 0.1;
end
if nargin < 8 || isempty(target_f0)
    target_f0 = [];
end
if nargin < 7 || isempty(harmonic_deviation)
    harmonic_deviation = 0.03;
end
if nargin < 6 || isempty(threshold_ratio)
    threshold_ratio = 0.7;
end
if nargin < 5 || isempty(threshold_amp)
    threshold_amp = 0.3;
end
if nargin < 4 || isempty(Fmax)
    Fmax = 800;
end
if nargin < 3 || isempty(Fmin)
    Fmin = 40;
end

f0 = 0;
acorr = [];
Nshifts = [];
Tshifts = [];

% Calculate index of min and max shift
Nshift_min = floor(1/(Fmax/Fs));
Nshift_max = ceil(1/(Fmin/Fs));
if Nshift_max+Nshift_min > length(A) || Nshift_min >= Nshift_max
    disp('Error in find_f0_timedomain2.m: Chunk size must be larger!')
    return
end

% Pull out the chunk of the signal and calculate a scale factor
% The scale factor is the probable maximum value
Achunk = A(1:Nshift_max);
scale_factor = sum(abs(Achunk))*2;

% Calculate a vector of all shifts
Nshifts = Nshift_min:Nshift_max;
Tshifts = Nshifts / Fs;

% Shift and calculate a makeshift autocorrelation
acorr = zeros(1,length(Nshifts));
index = 1;
for Nshift = Nshifts
    
    % Break out if we are shifting the chunk beyond the end of the vector A
    if Nshift_max + Nshift > length(A)
        Nshifts = Nshifts(1:index-1);
        Tshifts = Tshifts(1:index-1);
        acorr = acorr(1:index-1);
        break
    end
    
    % Calculate makeshift autocorrelation
    acorr(index) = 1-sum(abs(Achunk - A([1:Nshift_max] + Nshift)))./scale_factor; 
    index = index + 1;
    
end

% Try to make the autocorrelation "level" and with the minimum at zero
acorr = acorr - polyval(polyfit(Nshifts,acorr,1),Nshifts);
acorr = acorr - min(acorr);

% Find a list of all of the peaks above the threshold
[maxX,maxY,Imax,maxX_fit,peakX,peakY] = peakfind(Tshifts,acorr,5,1/Fmax,[],'',false);
% Keep only those harmonics that are within a certain height of the largest
% peak and whose height is above the threshold amplitude

% Make a counter to keep track of which peaks remain after each criteria is
% applied below...
Ipeak = 1:length(maxX_fit);

% Keep those whose heights are above the threshold
Ipeak = Ipeak(maxY > threshold_ratio*max(maxY) & maxY > threshold_amp);
maxX_fit = maxX_fit(Ipeak);

% If only 1 or zero peaks remain, return
N = length(Ipeak);
% if N == 1
%     f0 = 1./maxX_fit;
%     return
% elseif N == 0
%     return
% end
if N < 2
    return
end
    
% Keep only those harmonics whose frequency is less than a certain
% deviation around a multiple of the fundamental
[maxX_fit,Ix] = sort(maxX_fit);
Ipeak = Ipeak(Ix);
possible_f0 = zeros(1,N);
for n = 1:N-1
    if any(abs(1 - maxX_fit(n+1:N)./(2*maxX_fit(n))) < harmonic_deviation);
        if isempty(target_f0)
            f0 = 1./maxX_fit(n);
            Ipeak = Ipeak(n);
            break
        else
            possible_f0(n) = 1./maxX_fit(n);
        end
    end
end

% Now keep only the harmonic that is closest to the desired harmonic
if isempty(target_f0)
%    f0 = 1./maxX_fit(end);
%    return
else
    f0_error = abs(possible_f0./target_f0 - 1);
    [what,where] = min(f0_error);
    if what < target_tol
        f0 = possible_f0(where);
        Ipeak = Ipeak(where);
    end
end

% Now "fine tune" the frequency around the peak using a parabolic fit
if f0 > 0
    p = polyfit(peakX(Ipeak,:),peakY(Ipeak,:),2);
    Xfit = -p(2)/2/p(1);
    f0 = 1/Xfit;
end



function [maxX,maxY,Imax,maxX_fit,peakX,peakY] = peakfind(X,Y,Nmax,Xsep,Ymin,sortmethod,peakfit)


% This function finds the maxima of the given function, and
% returns the first Nvalues, sorted in decending order
%
% INPUTS:
%
% X = a vector of X-values for the signal
%
% Y = a fector of Y-values for the signal
%
% Nmax = maximum number of values to return
%
% Xsep = the minimum acceptable separation between "peaks"
% 
% Ymin = the minimum Y value to return
%
% sortmethod = how to sort the peaks...
%               'maxY' = Sort in decending order starting with maximum Y-valued peak
%               'minX' = Sort in ascending order starting with minimum X-valued peak
%
% peakfit = set to 'true' to do a parabolic fit and refine the maxX values
%
% OUTPUTS:
%
% maxX = the X-coordinates of all of the N peaks
%
% maxY = the Y-heights of all of the N peaks
%
% Imax = indices of all of the N peaks
%
% maxX_fit = a more refined list of the X-coordinates based on a parabolic
%           fit to each peak
%
% peakX = a Nx5 matrix with each row containing the 5 X-values around
%           each peak
%
% peakY = a Nx5 matrix with each row containing the 5 Y-values around each
%           peak

if nargin < 7 || isempty(peakfit)
    peakfit = true;
end
if nargin < 6 || isempty(sortmethod)
    sortmethod = 'maxY';
end
if nargin < 5
    Ymin = [];
end
if nargin < 4 || isempty(Xsep)
    Xsep = 0;
end
if nargin < 3 || isempty(Nmax)
    Nmax = 10;
end


% Differentiate, and find change in sign
Ydiff = diff(Y);
Imax1 = find(Ydiff(1:end-1) > 0 & Ydiff(2:end) < 0);
Imax2 = find(Ydiff(1:end-2) > 0 & Ydiff(2:end-1) == 0 & Ydiff(3:end) < 0);
Imax3 = find(Ydiff(1:end-3) > 0 & Ydiff(2:end-2) == 0 & Ydiff(3:end-1) == 0 & Ydiff(4:end) < 0);

% Concatenate all of the possible peaks
maxX = [X(Imax1+1) X(Imax2+1) X(Imax3+2)];
maxY = [Y(Imax1+1) Y(Imax2+1) Y(Imax3+2)];
Imax = [(Imax1+1) (Imax2+1) (Imax3+2)];

% Re-sort if desired
if strcmp(sortmethod,'minX')
    [maxX,IX] = sort(maxX,'ascend');
    maxY = maxY(IX);
    Imax = Imax(IX);
elseif strcmp(sortmethod,'maxY')
    [maxY,IY] = sort(maxY,'descend');
    maxX = maxX(IY);
    Imax = Imax(IY);
end

% Remove any peaks that are below the min peak value allowed
if ~isempty(Ymin)
    Irem = find(maxY < Ymin);
    maxY(Irem) = [];
    maxX(Irem) = [];
    Imax(Irem) = [];
end

% Remove any peaks that are closer than the minimum allowed peak seperation
% in X
if Xsep > 0
    Iremove = zeros(1,length(maxX));
    I = 0;
    for n = 2:length(maxX)
        if abs(maxX(n) - maxX(n-1)) < Xsep
            maxX(n) = maxX(n-1);
            maxY(n) = maxY(n-1);
            I = I+1;
            Iremove(I) = n;
        end
    end
    Iremove = Iremove(1:I);
    maxX(Iremove) = [];
    maxY(Iremove) = [];
    Imax(Iremove) = [];
end

% Remove any peaks beyond the max number to return
if length(maxX) > Nmax
    maxX = maxX(1:Nmax);
    maxY = maxY(1:Nmax);
    Imax = Imax(1:Nmax);
end

maxX_fit = maxX;
peakX = zeros(length(Imax),5);
peakY = zeros(length(Imax),5);
% Do a 2nd-order poly fit around the peak to pinpoint the peak location
for n = 1:length(Imax);
    if Imax(n) > 3 && Imax(n) < length(X)-2
        I = Imax(n)+[-2:2];
        if peakfit
            p = polyfit(X(I),Y(I),2);
            maxX_fit(n) = -p(2)/2/p(1);
        end
        peakX(n,1:5) = X(I);
        peakY(n,1:5) = Y(I);
    end
end

function [indices,freqs,notes,fund_index] = get_scale(scale)

indices = [];
freqs = [];
notes = {};
fund_index = [];

major_indices = [0 2 4 5 7 9 11];
minor_indices = [0 2 3 5 7 8 10];

fund_indices = [0 2 3 5 7 8 10];
fund_notes = {'A' 'B' 'C' 'D' 'E' 'F' 'G'};
scale_types = {'major','minor'};

if nargin < 1 || isempty(scale)
    % Return all possiblities of scales
    indices = {};
    for n = 1:length(fund_notes)
        for m = 1:2
            indices{end+1} = [fund_notes{n} scale_types{m}];
        end
    end
    return
else
    fund_index = fund_indices(strcmp(scale(1),fund_notes));
    if strcmp(scale(2:end),'major')
        indices = major_indices+fund_index;
    elseif strcmp(scale(2:end),'minor')
        indices = minor_indices+fund_index;
    else
        disp('ERROR in get_scale.m: unknown scale!')
        return
    end
end


indices = [indices-12 indices indices+12 indices+24 indices+36];
indices = indices(indices >=0 & indices <=48);

[indices,freqs,notes] = get_note_matrix(indices);

function [indices,freqs,notes] = get_note_matrix(note)

if nargin < 1
    note = '';
end

indices = (0:48)';
freqs = 55.*2.^(indices./12);
notes = {
    'A1'
    'A#/Bb1'
    'B1'
    'C1'
    'C#/Db1'
    'D1'
    'D#/Eb1'
    'E1'
    'F1'
    'F#/Gb1'
    'G1'
    'G#/Ab1'
    'A2'
    'A#/Bb2'
    'B2'
    'C2'
    'C#/Db2'
    'D2'
    'D#/Eb2'
    'E2'
    'F2'
    'F#/Gb2'
    'G2'
    'G#/Ab2'
    'A3'
    'A#/Bb3'
    'B3'
    'C3'
    'C#/Db3'
    'D3'
    'D#/Eb3'
    'E3'
    'F3'
    'F#/Gb3'
    'G3'
    'G#/Ab3'
    'A4'
    'A#/Bb4'
    'B4'
    'C4'
    'C#/Db4'
    'D4'
    'D#/Eb4'
    'E4'
    'F4'
    'F#/Gb4'
    'G4'
    'G#/Ab4'
    'A5'
    };

if ~isempty(note)
    if ischar(note)
        I = find(strcmp(notes,note));
        if ~isempty(I);
            indices = indices(I);
            freqs = freqs(I);
            notes = notes(I);
        else
            indices = [];
            freqs = [];
            notes = '';
        end
    else
        note = note(note>=0 & note<=48);
        indices = note;
        freqs = freqs(note+1);
        notes = notes(note+1);
    end
end





