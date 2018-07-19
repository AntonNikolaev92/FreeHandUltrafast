% External Client to receive AmfiTrack positions %

close all
clear all

addpath(genpath('G:\libs\matlab'));

port = 5555;
ipaddr = '10.47.9.240';
%ipaddr = '192.168.0.101';
szUSB = 64;
jTcpObj = [];

nFrames = [];
rcvBuf = [];

try
    jTcpObj = jtcp('request',ipaddr,port,'timeout',10000,'serialize',false);
    
    % receive position buffer
    while length(rcvBuf) < 4
        mssg = jtcp('read',jTcpObj);
        rcvBuf = [ rcvBuf, mssg ];
    end
    nFrames = typecast(rcvBuf(1:4),'uint32');
    nBytesToRead = nFrames*(4+szUSB)+4;
    fprintf('Frames to receive: %i (%i)\n', nFrames, nBytesToRead);
    while length(rcvBuf) < nBytesToRead
        mssg = jtcp('read',jTcpObj);
        rcvBuf = [ rcvBuf, mssg ];
    end
catch
    if ~isempty(jTcpObj), jtcp('close', jTcpObj); end;
    disp('something went wrong');
    return;
end

nStartByte = 5;
posBuf = rcvBuf(nStartByte: (nFrames*szUSB+nStartByte));
nStartByte = nStartByte+nFrames*szUSB;
posHex = dec2hex(typecast(posBuf,'uint8'));
tsBuf = rcvBuf(nStartByte:(length(rcvBuf)));

[ t, q ] = usbdata2pos(typecast(posBuf, 'uint8'), szUSB);
dt = typecast(tsBuf,'int32');

% make a time stamp
ts  = zeros(1,nFrames);
for iFrame = 2:nFrames
    ts(iFrame) = ts(iFrame-1)+double(dt(iFrame))*1e-9;
end

% Interpolate Coordinates with respect to acquisition frames
indexes = [];

for iFrame = 1:nFrames
    if ( t(1,iFrame) < -1e+4 ) || ( t(1,iFrame ) > 1e+4 )||...
          ( t(2,iFrame) < -1e+4 ) || ( t(2,iFrame ) > 1e+4 )|| ...
          ( t(3,iFrame) < -1e+4 ) || ( t(3,iFrame ) > 1e+4 )|| ...
          ( q(1,iFrame) < -1e+4 ) || ( q(1,iFrame ) > 1e+4 )|| ...
          ( q(2,iFrame) < -1e+4 ) || ( q(2,iFrame ) > 1e+4 )|| ...
          ( q(3,iFrame) < -1e+4 ) || ( q(3,iFrame ) > 1e+4 )
      indexes = indexes;
    else
        indexes = [indexes, iFrame];
    end
end

tmpt = t;
tmpq = q;
tmpts = ts;
t = ones(4,length(indexes));
q = zeros(4,length(indexes));
t(1,:) = tmpt(1,indexes);
t(2,:) = tmpt(2,indexes);
t(3,:) = tmpt(3,indexes);
t(4,:) = tmpt(4,indexes);
q(1,:) = tmpq(1,indexes);
q(2,:) = tmpq(2,indexes);
q(3,:) = tmpq(3,indexes);
q(4,:) = tmpq(4,indexes);
ts = tmpts(indexes);

subplot(3,3,1); plot(ts, t(1,:)); ylabel('tx'); grid on;
subplot(3,3,2); plot(ts, t(2,:)); ylabel('ty'); grid on;
subplot(3,3,3); plot(ts, t(3,:)); ylabel('tz'); grid on;
subplot(3,3,4); plot(ts, q(1,:)); ylabel('qx'); grid on;
subplot(3,3,5); plot(ts, q(2,:)); ylabel('qy'); grid on;
subplot(3,3,6); plot(ts, q(3,:)); ylabel('qz'); grid on;
subplot(3,3,7); plot(ts); grid on;

