function [rf_d] = rfDepcon(rf, dt, rayParam, dz, model, normalize)
% RFDEPCON Depth convert receiver functions using a one-dimensional
%          earth model.
%
% >> [rf_d] = RFDEPCON(rf, dt, rayParameter, dz, model)
%
%---Input Variables------------------------------------------------
% rf        - Receiver function data
% dt        - Sample rate of data (s)
% rayParam  - P-wave ray parameter (s/km)
% dz        - Depth sampling rate (km)
% model     - Earth model (accepted values: 'iasp91', 'prem')
% normalize - Normalize the receiver function amplitude ('true','false')
%
%---Output Variables-----------------------------------------------
% rf_d      - Depth converted receiver function
%
%------------------------------------------------------------------
% Last updated 6/02/2020 by aburky@princeton.edu
%------------------------------------------------------------------

% Check the input Earth model
if strcmp(model, 'iasp91')
    emodel = dlmread('/Users/aburky/IFILES/MODELS/IASP91/IASP91.tvel','',2,0);
    z = emodel(1:125,1);
    vp = emodel(1:125,2);
    vs = emodel(1:125,3);
elseif strcmp(model, 'prem')
    emodel = dlmread('/Users/aburky/IFILES/MODELS/PREM_ISO.tvel','',2,0);
    z = emodel(1:42,1);
    vp = emodel(1:42,2);
    vs = emodel(1:42,3);
else
    error(['Invalid Earth model. Supported values are '...
           '''iasp91'' and ''prem''']);
end

% Check the input ray parameter
if rayParam > 1
    error('Invalid ray parameter, units should be (s/km)')
else
    P = rayParam;
end

% Find discontinuities in the Earth model
idx = find(diff(z)==0) + 1;

% Handle discontinuities in the Earth model
for i=1:length(idx)+1
    if i == 1
        zp{i} = z(i:idx(i)-1);
        vpp{i} = vp(i:idx(i)-1);
        vsp{i} = vs(i:idx(i)-1);
    elseif i > 1 && i <= length(idx)
        zp{i} = z(idx(i-1):idx(i)-1);
        vpp{i} = vp(idx(i-1):idx(i)-1);
        vsp{i} = vs(idx(i-1):idx(i)-1);
    elseif i == length(idx)+1
        zp{i} = z(idx(i-1):end);
        vpp{i} = vp(idx(i-1):end);
        vsp{i} = vs(idx(i-1):end);
    end
end

% Construct interpolants
for i = 1:length(idx)
    vp_int{i} = interp1(zp{i},vpp{i},zp{i}(1):dz:zp{i}(end));
    vs_int{i} = interp1(zp{i},vsp{i},zp{i}(1):dz:zp{i}(end));
    zzp{i} = zp{i}:dz:zp{i}(end);
end

% Construct depth conversion integrand terms (Chevrot et al.)
for j = 1:length(P)
    for i = 1:length(idx)
        a{i}{j} = sqrt((vs_int{i}.^(-2))-P(j)^2*...
                      ((6371-zzp{i})./6371).^(-2));
        b{i}{j} = sqrt((vp_int{i}.^(-2))-P(j)^2*...
                      ((6371-zzp{i})./6371).^(-2));
        v{i}{j} = a{i}{j} - b{i}{j};
    end
end

% Solve the integral piecewise
for j = 1:length(P)
    for i = 1:length(idx)
        tz{j}{i}(1) = 0;
        for k = 2:length(zzp{i})
            tz{j}{i}(k) = trapz(zzp{i}(1:k),v{i}{j}(1:k));
        end
    end
end

% Concatenate into a time to depth conversion vector
for j = 1:length(P)
    ttd{j} = tz{j};
    for i = 1:length(idx)-1
        ttd{j}{i+1} = ttd{j}{i+1} + ttd{j}{i}(end);
    end
end

for j = 1:length(P)
    ttz{j} = {[ ttd{j}{:} ]};
    ttzm(:,j) = cell2mat(ttz{j});
end

% Delete duplicate entires in ttzm (corresponding to discontinuities)
dup_idx = find(diff(ttzm(:,1)) == 0);
ttzm(dup_idx,:) = [];

% Keep real components of ttzm
if isreal(ttzm) == 0
    ttzm = real(ttzm);
    warning('Depth conversion matrix has imaginary values')
end

% Normalize receiver function amplitude by max in small time window
if strcmp(normalize,'true')
    rf = rf/max(rf(1:fix(12/dt)));
    % rf = rf/max(abs(rf));
    % Align maximum to zero time
    [~,shidx] = max(rf(1:fix(12/dt)));
    % [~,shidx] = max(abs(rf));
else
    % Default time shift of 10 seconds
    shidx = fix(10/dt);
end

% Remove 10 second time shift
% shidx = round(10/dt);
rf = rf(shidx:end);

% Depth convert the data
ttzm_r = round(ttzm/dt);
% Shift all data by one, to deal with ttzm(1,:) being equal to 0
ttzm_r(1,:) = 1;

rf_d = rf(ttzm_r);


