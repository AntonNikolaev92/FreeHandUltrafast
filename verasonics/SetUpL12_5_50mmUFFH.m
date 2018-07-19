

clear all

filePath = 'G:\murab_ugb\Free Hand Scanning\FreehandUltrafast\Verasonics';

nPackages = 4;  % number of packeges to transfer to host
nFrames = 1000; % number frames in package

P.startDepth = 2;   % Acquisition depth in wavelengths
P.endDepth = 192;   % This should preferrably be a multiple of 128 samples.

% Define system parameters.
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;

% Specify Trans structure array.
Trans.name = 'L12-5 50mm';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % L12-5_50mm transducer is 'known' transducer so we can use computeTrans.

% Specify PData structure array.
PData.PDelta = [Trans.spacing, 0, 0.5];
PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3)); % startDepth, endDepth and pdelta set PData.Size.
PData.Size(2) = ceil((128*Trans.spacing)/PData.PDelta(1));
PData.Size(3) = 1;      % single image page
PData.Origin = [-Trans.spacing*(128-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096; % Rows per package (super frame)
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 40;        % 40 frames used for RF cineloop.
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 20;
Resource.DisplayWindow(1).Title = 'L12-5_50mmFlashs';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
    DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = 4096*nFrames; % Rows per package (super frame)
Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(2).numFrames = nPackages;        % 40 frames used for RF cineloop.

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,1,1];   % A, B, C, D  for 8.1 MHz transmit frequency

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
    'Origin', [0.0,0.0,0.0], ...
    'aperture', 65, ...
    'Apod', ones(1,Resource.Parameters.numTransmit), ...
    'focus', 0.0, ...
    'Steer', [0.0,0.0], ...
    'Delay', zeros(1,Resource.Parameters.numTransmit)), 1, 1);

% Specify TGC Waveform structure.
TGC.CntrlPts = [234,368,514,609,750,856,1023,1023];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays -
%   endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.
%   InputFilter - The same coefficients are used for all channels. The
%              coefficients below give a broad bandwidth bandpass filter.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', zeros(1,128), ...
    'aperture', 65, ...
    'startDepth', P.startDepth, ...
    'endDepth',maxAcqLength, ...
    'TGC', 1, ...
    'bufnum', 1, ...
    'framenum', 1, ...
    'acqNum', 1, ...
    'sampleMode', 'NS200BW', ...
    'mode', 0, ...
    'callMediaFunc', 0),1,nPackages*nFrames + Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes 
for iPackage = 1:nPackages  % for UF acquisition
    i = (iPackage-1)*nFrames;
    for iFrame = 1:nFrames
        Receive(i+iFrame).callMediaFunc = 1;
        Receive(i+iFrame).Apod(1:128) = 1.0;
        Receive(i+iFrame).aperture = 65;
        Receive(i+iFrame).framenum = iPackage;
        Receive(i+iFrame).acqNum = iFrame;
    end
end
for iFrame = 1:Resource.RcvBuffer(1).numFrames % for normal acquisition
    Receive(nPackages*nFrames + iFrame)
    Receive(i+iFrame).callMediaFunc = 1;
    Receive(i+iFrame).Apod(1:128) = 1.0;
    Receive(i+iFrame).aperture = 65;
    Receive(i+iFrame).framenum = iFrame;
    Receive(i+iFrame).acqNum = 1;
end

% Specify Recon structure arrays.
Recon = struct('senscutoff', 0.6, ...
    'pdatanum', 1, ...
    'rcvBufFrame', -1, ...     % use most recently transferred frame
    'IntBufDest', [1,1], ...
    'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
    'RINums', 1);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 0, ...  % replace IQ data.
    'txnum', 1, ...
    'rcvnum', 1, ...
    'scaleFactor', 2.0, ...
    'regionnum', 1), 1, 1);

% Specify Process structure array.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
    'framenum',-1,...   % (-1 => lastFrame)
    'pdatanum',1,...    % number of PData structure to use
    'pgain',1.0,...            % pgain is image processing gain
    'reject',2,...      % reject level
    'persistMethod','simple',...
    'persistLevel',pers,...
    'interpMethod','4pt',...  %method of interp. (1=4pt)
    'grainRemoval','none',...
    'processMethod','none',...
    'averageMethod','none',...
    'compressMethod','power',...
    'compressFactor',40,...
    'mappingMethod','full',...
    'display',1,...      % display image after processing
    'displayWindow',1};
                       
Process(2).classname = 'External';
Process(2).method = 'extSimple';
Process(2).Parameters = {'srcbuffer','none',...  % name of buffer to process.
                         'srcbufnum',0,...
                         'srcframenum',0,...
                         'dstbuffer','none'};

%% Describe UFFH

% Specify SeqControl structure arrays.
scJump =1;
SeqControl(scJump).command = 'jump'; % jump back to start.
SeqControl(scJump).argument = 1;
scTimeToNextAcq = 2;
SeqControl(scTimeToNextAcq).command = 'timeToNextAcq'; % time between synthetic aperture acquisitions
SeqControl(scTimeToNextAcq).argument = 200; % 200 usec (5 kfps)
scTrigOut = 3;
SeqControl(scTrigOut).command = 'triggerOut';
scLoopCnt0 = 4;
SeqControl(scLoopCnt0).command = 'loopCnt';
SeqControl(scLoopCnt0).argument = 0;        
scLoopCnt1 = 5;
SeqControl(scLoopCnt1).command = 'loopCnt';
SeqControl(scLoopCnt1).argument = 1;
scLoopTst = 6;
SeqControl(scLoopTst).command = 'loopTst';
scLoopTst1 = 7;
SeqControl(scJumpToBegin).command = 'jump';
SeqControl(scJumpToBegin).argument = 1;

nsc = 8; % nsc is count of SeqControl objects

n = 1; % n is count of Events

Event(n) = createEvent('initialise case 0');
Event(n).seqControl = scLoopCnt0;
n = n+1;
Event(n) = createEvent('jump to start');
Event(n).seqControl = scJump;
n = n+1;
Event(n) = createEvent('initialise case 1');
Event(n).seqControl = scLoopCnt1;
n = n+1;

% switch:
Event(n) = createEvent('check fist condition');
Event(n).seqControl = scLoopTst0;
SeqControl(scLoopTst0).argument = n+3;
n = n+1;

% case 0:
for i = 1:Resource.RcvBuffer(1).numFrames       % Acquire frame
    Event(n).info = 'acquire ray line.';
    Event(n).tx = 1;   % use next TX structure.
    Event(n).rcv = 1;   
    Event(n).recon = 0;      % no reconstruction.
    Event(n).process = 0;    % no processing
    Event(n).seqControl = 2;
    n = n+1;
    Event(n-1).seqControl = [3,nsc]; % modify last acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recon and process'; 
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 1;      % reconstruction
    Event(n).process = 1;    % process
    Event(n).seqControl = 0;
    if (floor(i/5) == i/5)&&(i ~= Resource.RcvBuffer(1).numFrames) % Exit to Matlab every 5th frame 
        Event(n).seqControl = 4; % return to Matlab
    end
    n = n+1;
end
Event(n) = createEvent('restore loop control');
Event(n).seqControl = scLoopCnt0;
n = n+1;
Event(n) = createEvent('jump back to start');
Event(n).seqControl = scJump;
n = n+1;

%case 1:
Event(n) = createEvent('start UF acquisition trigger');
Event(n).seqControl = scTrigOut;
n = n+1;
for iPackage = 1:nPackages
    % Collect a singel package to Transfer
    i = (iPackage-1)*nFrames;
    for iFrame = 1:nFrames
        Event(n).info = 'aperture';
        Event(n).tx = 1;       % transmit from the single aperture
        Event(n).rcv = i+iFrame;      % receive with the single aperture
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        if iFrame == nFrames
            Event(n).seqControl = [scTimeToNextAcq, nsc];
            SeqControl(nsc).command = 'transferToHost';
            nsc = nsc + 1;
        else
            Event(n).seqControl = [scTimeToNextAcq] ; % time between frames, SeqControl struct defined below.
        end
        n = n+1;
    end
end
%case 1:
Event(n) = createEvent('stop UF acquisition trigger');
Event(n).seqControl = scTrigOut;
n = n+1;
Event(n) = createEvent('call external function');
Event(n).process = 2;
n = n+1;
Event(n) = createEvent('restore loop control');
Event(n).seqControl = scLoopCnt1;
n = n+1;
Event(n) = createEvent('jump back to start');
Event(n).seqControl = scJumpToBegin;


% User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
    'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
    'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutoffCallback');

% - Range Change
MinMaxVal = [64,300,P.endDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
    end
end
UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
    'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%RangeChangeCallback');
UI(3).Control = {'UserC2','Style','VsPushButton','Label', 'Acquire'};
UI(3).Callback = text2cell('%AcquirePositionCallback');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 2;

% Save all the structures to a .mat file.
filename = [ filePath, '/', 'MatFiles', '/', 'L12-5_50mmUFFH', '.mat'];
save(filename);
disp('Data stored. Start Acting');
VSX
return

% **** Callback routines to be converted by text2cell function. ****

%AcquirePositionCallback
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Parameters',1,'startEvent',3};
assignin('base','Control', Control);
return
%AcquirePositionCallback

%SensCutoffCallback - Sensitivity cutoff change
ReconL = evalin('base', 'Recon');
for i = 1:size(ReconL,2)
    ReconL(i).senscutoff = UIValue;
end
assignin('base','Recon',ReconL);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'Recon'};
assignin('base','Control', Control);
return
%SensCutoffCallback

%RangeChangeCallback - Range change
simMode = evalin('base','Resource.Parameters.simulateMode');
% No range change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.endDepth'));
    return
end
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

P = evalin('base','P');
P.endDepth = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        P.endDepth = UIValue*scaleToWvl;
    end
end
assignin('base','P',P);

evalin('base','PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));');
evalin('base','PData(1).Region = computeRegions(PData(1));');
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
Receive = evalin('base', 'Receive');
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
for i = 1:size(Receive,2)
    Receive(i).endDepth = maxAcqLength;
end
assignin('base','Receive',Receive);
evalin('base','TGC.rangeMax = P.endDepth;');
evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
evalin('base','if VDAS==1, Result = loadTgcWaveform(1); end');
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback
