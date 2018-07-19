function newEvent = createEvent( info )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    newEvent.info = info;
    newEvent.tx = 0;   % use next TX structure.
    newEvent.rcv =0;   
    newEvent.recon = 0;      % no reconstruction.
    newEvent.process = 0;    % no processing
    newEvent.seqControl = 0;
end

