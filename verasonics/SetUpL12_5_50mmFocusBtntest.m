% Notice: 
% File name
% Description: 
% Last update 

clear all

P.startDepth = 0;  % Define P.startDepth and P.endDepth at top for use in defining other parameters.
P.endDepth = 208;
P.numRays = 128;
P.numTx = 24;    % number of active transmitters in TX aperture.
P.txFocus = 128;   % transmit focal pt in wavelengths

% Define system parameters.
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
Resource.Parameters.startEvent = 1;

% Specify Trans structure array.
Trans.name = 'L12-5 50mm';
Trans.units = 'wavelengths'; % required in Gen3 to prevent default to mm units
% note nominal center frequency in computeTrans is 7.813 MHz
Trans.maxHighVoltage = 35;  % set maximum high voltage limit for pulser supply.
Trans = computeTrans(Trans);  % L28-15 transducer is 'known' transducer so we can use computeTrans.

% Specify PData structure array.
PData.PDelta = [Trans.spacing, 0, 0.5];
PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3)); % startDepth, endDepth and pdelta set PData.Size.
PData.Size(2) = ceil((Trans.numelements*Trans.spacing)/PData.PDelta(1));
PData.Size(3) = 1;      % single image page
PData.Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% Define numRays rectangular regions centered on TX beam origins.
PData(1).Region = repmat(struct('Shape',struct('Name','Rectangle',...
                                               'Position',[0,0,P.startDepth],...
                                               'width',Trans.spacing,...
                                               'height',P.endDepth-P.startDepth)),1,P.numRays);
firstRayLocX = -((P.numRays-1)/2)*Trans.spacing;
for i = 1:P.numRays
    PData.Region(i).Shape.Position(1) = firstRayLocX + (i-1)*Trans.spacing;
end
PData.Region = computeRegions(PData);

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = P.numRays*208*4*2; % this size allows for all rays, with range of up to 400 wvlngths
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 20;        % 20 frames used for RF cineloop.
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L12-5_50mm_256RyLns';
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

% Specify Transmit waveform structure.  
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,2,1];

% Specify TX structure array.  
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'aperture', 1, ...
                   'Apod', zeros(1,Resource.Parameters.numTransmit), ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Resource.Parameters.numTransmit)), 1, P.numRays);

scaleToWvl = 1;
if strcmp(Trans.units, 'mm')
    scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
end              
               
% - Set event specific TX attributes.
for n = 1:P.numRays  % specify P.numRays transmit events
    % Set transmit Origins.
    TX(n).Origin = [(firstRayLocX + (n-1)*Trans.spacing), 0.0, 0.0];
    % Compute available transmit mux aperture
    [Dummy,ce] = min(abs(scaleToWvl*Trans.ElementPos(:,1)-TX(n).Origin(1))); % ce is closest ele to cntr of aper.
    % Determine first element of mux aperture.
    lft = round(ce - 64);
    if lft < 1, lft = 1; end
    if lft > 129, lft = 129; end
    TX(n).aperture = lft;
    % Compute TX.Apod within mux aperture.
    lft = round(ce - P.numTx/2);
    if lft < 1, lft = 1; end
    rt = round(ce + P.numTx/2);
    if rt > 256, rt = 256; end
    TX(n).Apod((lft-(TX(n).aperture-1)):(rt-(TX(n).aperture-1))) = 1;
    % Apply apodization function.
    [RIndices,CIndices,V] = find(TX(n).Apod);
    V = kaiser(size(V,2),1);
    TX(n).Apod(CIndices) = V;
    % Compute transmit delays
    TX(n).Delay = computeTXDelays(TX(n));
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [395,535,650,710,770,830,890,950];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays 
RcvProfile.LnaZinSel = 31;
%   P.endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
maxAcqLength = P.endDepth;
Receive = repmat(struct('aperture', 1, ...
                        'Apod', zeros(1,128), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'BS100BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0),1,P.numRays*Resource.RcvBuffer(1).numFrames);
                    
% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    Receive(P.numRays*(i-1)+1).callMediaFunc = 1; 
    for j = 1:P.numRays
        Receive(P.numRays*(i-1)+j).aperture = TX(j).aperture; % mux aperture same as transmit
        Receive(P.numRays*(i-1)+j).Apod(1:128) = 1.0;
        Receive(P.numRays*(i-1)+j).framenum = i;
        Receive(P.numRays*(i-1)+j).acqNum = j;
    end    
end

% Specify Recon structure arrays.
% - We need one Recon structure which will be used for each frame. 
Recon = struct('senscutoff', 0.5, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'ImgBufDest', [1,-1], ...
               'RINums',1:P.numRays);

% Define ReconInfo structures.
% We need 2*na ReconInfo structures for 2-1 synthetic aperture and na steering angles.
ReconInfo = repmat(struct('mode', 0, ...  % default is to accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 0), 1, P.numRays);
% - Set specific ReconInfo attributes.
for j = 1:P.numRays
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = j;
end

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
Process(2).method = 'myProc';
Process(2).Parameters = {'srcbuffer','none',...  % name of buffer to process.
  'srcbufnum',0,...
  'srcframenum',0,...
  'dstbuffer','none'};

                     
% Specify SeqControl structure arrays.
%  - Jump back to start.
scJump = 1;
SeqControl(scJump).command = 'jump';
SeqControl(scJump).argument = 4;
scTimeToNextAcq = 2;
SeqControl(scTimeToNextAcq).command = 'timeToNextAcq';  % time between ray line acquisitions
SeqControl(scTimeToNextAcq).argument = 400; % 400 usec
scTimeToNextAcqFrame = 3;
SeqControl(scTimeToNextAcqFrame).command = 'timeToNextAcq';
SeqControl(scTimeToNextAcqFrame).argument = 10000;   % 10000 usec = 10msec time between frames
scReturnToMatlab = 4;
SeqControl(scReturnToMatlab).command = 'returnToMatlab';
scLoopCnt0 = 5;
SeqControl(scLoopCnt0).command = 'loopCnt';
SeqControl(scLoopCnt0).argument = 0;        
scLoopCnt1 = 6;
SeqControl(scLoopCnt1).command = 'loopCnt';
SeqControl(scLoopCnt1).argument = 1;
scLoopTst0 = 7;
SeqControl(scLoopTst0).command = 'loopTst';
scLoopTst1 = 8;
SeqControl(scLoopTst1).command = 'loopTst';
scJumpToBegin = 9;
SeqControl(scJumpToBegin).command = 'jump';
SeqControl(scJumpToBegin).argument = 1;

nsc = 10; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;

Event(n) = createEvent('initialise case 0');
Event(n).seqControl = scLoopCnt0;
n = n+1;
Event(n) = createEvent('jump to start');
Event(n).seqControl = scJump;
n = n+1;
Event(n) = createEvent('initialise case 1');
Event(n).seqControl = scLoopCnt1;
n = n+1;

%SeqControl(scJump).argument = n;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:P.numRays        % Acquire frame
        Event(n).info = 'acquire ray line.';
        Event(n).tx = j;   % use next TX structure.
        Event(n).rcv = P.numRays*(i-1)+j;   
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 2;
        n = n+1;
    end
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

% switch:
Event(n) = createEvent('check fist condition');
Event(n).seqControl = scLoopTst0;
SeqControl(scLoopTst0).argument = n+3;
n = n+1;

% case 0:
Event(n) = createEvent('restore loop control');
Event(n).seqControl = scLoopCnt0;
n = n+1;
Event(n) = createEvent('jump back to start');
Event(n).seqControl = scJump;
n = n+1;

%case 1:
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
UI(1).Callback = text2cell('%SensCutOffCallback');

% - Range Change
wls2mm = 1;
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        AxesUnit = 'mm';
        wls2mm = Resource.Parameters.speedOfSound/1000/Trans.frequency;
    end
end
UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',[64,300,P.endDepth]*wls2mm,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%RangeChangeCallback');
             
% - Transmit focus change
UI(3).Control = {'UserB4','Style','VsSlider','Label',['TX Focus (',AxesUnit,')'],...
                 'SliderMinMaxVal',[20,320,P.txFocus]*wls2mm,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(3).Callback = text2cell('%txFocusCallback');
             
% - TX Aperture change
UI(4).Control = {'UserB3','Style','VsSlider','Label','TX Aper',...
                 'SliderMinMaxVal',[1,128,P.numTx],'SliderStep',[1/128,1/64],'ValueFormat','%3.0f'};
UI(4).Callback = text2cell('%ApertureCallback');

UI(5).Control = {'UserC2','Style','VsPushButton','Label', 'TakeFrame'};
UI(5).Callback = text2cell('%AcquirePositionCallback');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 5;

% Save all the structures to a .mat file.
filename = 'MatFiles/L12-5_50mm_128Foc';
save(filename);
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

%SensCutOffCallback - Sensitivity cutoff change
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
%SensCutOffCallback

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

PData = evalin('base','PData');
PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3)); % P.startDepth, P.endDepth and pdelta set PData.Size(1).
for i = 1:size(PData.Region,2)
    PData.Region(i).Shape.height = P.endDepth-P.startDepth;
end
PData.Region = computeRegions(PData);
assignin('base','PData',PData);
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData.Size(1)*PData.PDelta(3)/Resource.DisplayWindow(1).pdelta);');
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

%txFocusCallback
simMode = evalin('base','Resource.Parameters.simulateMode');
% No focus change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.txFocus'));
    return
end
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

P = evalin('base','P');
P.txFocus = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        P.txFocus = UIValue*scaleToWvl;    
    end
end
assignin('base','P',P);

% - Redefine event specific TX attributes for the new focus.
TX = evalin('base', 'TX');
for n = 1:P.numRays
    % write new focus value to TX
    TX(n).focus = P.txFocus;
    TX(n).Delay = computeTXDelays(TX(n));
end
assignin('base','TX', TX);
% Set Control command to update TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);
return
%txFocusCallback

%ApertureCallback
simMode = evalin('base','Resource.Parameters.simulateMode');
% No focus change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.numTx'));
    return
end
Trans = evalin('base', 'Trans');
P = evalin('base','P');
P.numTx = UIValue;
assignin('base','P',P);

TX = evalin('base', 'TX');
scaleToWvl = evalin('base','scaleToWvl');
rayDelta = Trans.numelements*Trans.spacing/P.numRays;
firstRayLocX = -((Trans.numelements-1)/2)*Trans.spacing;
% - Redefine event specific TX attributes for the new aperture.
for n = 1:P.numRays  % specify P.numRays transmit events
    % Set transmit Origins.
    TX(n).Origin = [(firstRayLocX + (n-1)*rayDelta), 0.0, 0.0];
    % Compute available transmit mux aperture
    [Dummy,ce] = min(abs(scaleToWvl*Trans.ElementPos(:,1)-TX(n).Origin(1))); % ce is closest ele to cntr of aper.
    lft = round(ce - 64);
    if lft < 1, lft = 1; end
    if lft > 129, lft = 129; end
    TX(n).aperture = lft;
    % Compute TX.Apod within mux aperture.
    lft = round(ce - P.numTx/2);
    if lft < 1, lft = 1; end
    rt = round(ce + P.numTx/2);
    if rt > 256, rt = 256; end
    TX(n).Apod((lft-(TX(n).aperture-1)):(rt-(TX(n).aperture-1))) = 1;
    % Apply apodization function.
    [RIndices,CIndices,V] = find(TX(n).Apod);
    V = kaiser(size(V,2),1);
    TX(n).Apod(CIndices) = V;
    % Compute transmit delays
    TX(n).Delay = computeTXDelays(TX(n));
end
assignin('base','TX', TX);
% Set Control command to update TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);
return
%ApertureCallback