function saveRF(rf,directory)
% SAVERF Saves receiver function data to a SAC file and fills all of the 
%        appropriate header variables.
%
% >> saveRF(rf,directory)
%
%---Input Variables--------------------------------------------------------
% rf        - structure containing the receiver function data
% directory - string indicating the path to the folder where you would like
%             to save the data
%
%---External Dependencies--------------------------------------------------
% F_SAC: fwrite_sac
% qcRF      - function which calculates the receiver function quality
%             control parameter
%
%--------------------------------------------------------------------------
% Last updated 1/5/2021 by aburky@princeton.edu
%--------------------------------------------------------------------------

% Receiver function timeseries data
sacmat.data = rf.d;
sacmat.hdr = rf.h;

% Save the Vertical and Radial SNR to SAC Header Variables (USER0 + USER1)
sacmat.hdr.user(1) = rf.vsnr;
sacmat.hdr.user(2) = rf.rsnr;

% Save the RMS error to a SAC Header Variable (USER2)
sacmat.hdr.user(3) = 100 * (1 - rf.rms(end));

% Save the Quality Control parameter (nu) to a SAC Header Variable (USER3)
sacmat.hdr.user(4) = qcRF(rf.d,rf.t);

% Format the receiver function filename
[~, mm, dd, ~, ~, ~] = datevec(datenum(sacmat.hdr.nzyear,...
                               1,sacmat.hdr.nzjday));
if strcmp(sacmat.hdr.khole,'-12345')
    locID = '';
elseif strcmp(sacmat.hdr.khole,'')
    locID = '';
else
    locID = sacmat.hdr.khole;
end
fname = sprintf('%i.%02d.%02d.%02d.%02d.%02d.%s.%s.%s.RF.SAC',...
                sacmat.hdr.nzyear,mm,dd,sacmat.hdr.nzhour,...
                sacmat.hdr.nzmin,sacmat.hdr.nzsec,sacmat.hdr.knetwk,...
                sacmat.hdr.kstnm,locID);
            
% Save the receiver function data as a SAC file
fwrite_sac(sacmat,fullfile(directory,fname));
