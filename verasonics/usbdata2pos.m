function [ t, q ] = usbdata2pos( data, szPackage )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

nFrames = floor(length(data)/szPackage);

t = zeros(3,nFrames);
q = zeros(4,nFrames);

for iFrame = 1:nFrames
    iStart = (iFrame-1)*szPackage + 1;
    
    % XPos
    tx = amfiData2Dec(data(iStart+10:iStart+12)'); tx = tx*0.01;% mm
    ty = amfiData2Dec(data(iStart+13:iStart+15)'); ty = ty*0.01;% mm
    tz = amfiData2Dec(data(iStart+16:iStart+18)'); tz = tz*0.01;% mm
    t(:,iFrame) = [tx,ty,tz]';
    qx = amfiData2Dec(data(iStart+19:iStart+21)'); qx = qx*1e-6;
    qy = amfiData2Dec(data(iStart+22:iStart+24)'); qy = qy*1e-6;
    qz = amfiData2Dec(data(iStart+25:iStart+27)'); qz = qz*1e-6;
    qw = amfiData2Dec(data(iStart+28:iStart+30)'); qw = qw*1e-6;
    q(:,iFrame)  = [qw, qx, qy, qz ]';


end

end

function [ decNumber] = amfiData2Dec(values)
    lenData = length(values);
    if lenData~=3, decNumber = 0; return; end;
    % Format Number
    decNumber = double(values(3))*16.^4 + double(values(2))*16.^2 + double(values(1));
    numTh = pow2(23) - 1 ;
    numMax = pow2(24) - 1;
    
    if decNumber > numTh, decNumber = decNumber - numMax; end;
end
