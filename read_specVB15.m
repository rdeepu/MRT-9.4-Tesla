%% Spectroscopic Measurements VB15
%  Calculates either T1, T2 or the Reference Voltage, depending on name of
%  fname. If path is set to 1, it expects fname to be a path in which the
%  user should select the corresponding file.

function [sig,spec, fitVal, mrprot] = read_specVB15(fname, path)
%% 
%  INPUT:
%
% * fname: Name of input file or path (depending on path)
% * (path): If set to 1, fname is expected to be a path and not a name
%% Read Data and get Protocol
if nargin < 1 
    info = 'Please select binary file to read';
    [fname,pathname]=uigetfile('*.dat',info);  
    fname = [pathname fname];
elseif nargin < 3 && path == 1
    info = 'Please select binary file to read';
    [fname,pathname]=uigetfile('',info, fname);
    fname = [pathname fname];
end

os=2;
mrprot = readVB17Header(fname);
if  ~isstruct(mrprot)
    errordlg('Error while reading VB17Header. Not able to continue');
    sig = -1;
    spec = -1;
    fitVal = -1;
    return;
end
nCol = mrprot.MeasYaps.sSpecPara.lVectorSize*os;
nCha = mrprot.Meas.iMaxNoOfRxChannels;
nAve = mrprot.Meas.NAveMeas;
nMea = mrprot.Meas.NRepMeas; % removed +1

dw    = mrprot.MeasYaps.sRXSPEC.alDwellTime{1}*1e-9;
bw          = (mrprot.MeasYaps.sRXSPEC.alDwellTime{1}*1e-9)^-1;

if ~isempty(findstr(fname,'T1'))
    disp('T1 Measurement')
    mode = 'T1';
    first        = mrprot.MeasYaps.sWiPMemBlock.adFree{2};
    delta        = mrprot.MeasYaps.sWiPMemBlock.adFree{3};   
    img = nMea;
elseif ~isempty(findstr(fname,'T2'))
    disp('T2 Measurement')
    mode = 'T2';  
    delta        = mrprot.MeasYaps.sWiPMemBlock.adFree{1};
    first        = mrprot.Meas.alTE(1)/1000;
    img          = 1;
elseif ~isempty(findstr(fname,'FA')) || ~isempty(findstr(fname,'RefVolt'))
    disp('FA Calculation')
    mode = 'FA';
    first = [];
    delta = [];
    img = nMea/2;
else
    mode = '';
    disp('Not supported')
end


fid=fopen(fname);
fseek(fid,0,'bof');
firstLine=fread(fid,1,'uint32'); % first 4 Bytes tell where data section starts
fseek(fid,firstLine,'bof');

% allocate memory for data matrix
data=zeros(nCol,nAve,nMea,nCha);
voltage = zeros(1, nMea);

counter=0;
for k=1:nMea
    for m=1:nAve
        for l=1:nCha
            counter=counter+1;
            mdh = read_mdh(fid);
            if strcmp(mode,'FA')
                voltage(k) = mdh.sLoopCounter.ushPhase;
            end
            tmp = fread(fid,2*nCol,'float'); %factor 2, because of (real,imag)
            data(:,m,k,l)=tmp(1:2:2*nCol)+1i*tmp(2:2:2*nCol);  
            fseek(fid,firstLine+counter*(nCol*2*4+128),-1); % 4 Bytes real and imag (2)
        end
    end
end

%% Transform Data
% Apply iFFT, perform averaging, and calculate SOS if more than one recieve channel was used

sig = squeeze(sum(data,2));
spec=ifftshift(ifft(ifftshift(sig)));
spec=makesos(spec,3);

%% Plot Time Domain Signal and Spectrum

ts = dw *linspace(0,nCol,nCol);
fs = bw*linspace(-0.5,0.5,nCol);
figure,plot(ts,real(sig(:,img)),'k')
    title('Signal in Time Domain (Real Part)')
    xlabel('t [s]')
    ylabel('Signal [a.u.]')
    grid on

figure,plot(fs,abs(spec(:,6)),'k')
    title('Spectrum (Absolute Value)')
    xlabel('\nu [Hz]')
    ylabel('Signal [a.u.]')
    grid on

%% Fit Data
specMax=max(abs(spec));
time=first:delta:(first+delta*nMea-1);

if strcmp(mode,'T1')
    [value,inv_time]=min(specMax);
    est_T1=time(inv_time)./log((specMax(1)+specMax(nMea))./specMax(1));
    est_M0=specMax(end);
    est_Eff=1;
    xout=fminsearch('t1_3pabs_eff',[est_T1 est_M0 est_Eff],[],time,specMax);

    figure, plot(time,specMax,'k+');hold on;
    t=[1:0.1:max(time)].';plot(t,abs(xout(2)*(1-(1+xout(3))*exp(-t./xout(1)))),'k');
    title('T_1 Measurement')
    xlabel('TI [ms]')
    ylabel('Signal [a.u.]')
    legend('Measurement','Fit','Location','Best')
    grid on
    disp(['fitted T1 : ' num2str(xout(1)) ' ms']);
    disp(['fitted M0 : ' num2str(xout(2)*1e5 ) ' a.u.']);
    disp(['Efficiency: ' num2str(xout(3)) ]);
    fitVal = xout(1);
elseif strcmp(mode,'T2')
    est_T2=5;
    est_M0=specMax(1);
    est_C=specMax(end);

    xout=fminsearch('t2_3p',[est_T2 est_M0 est_C],[],time,specMax);

   figure, plot(time,specMax,'k+');hold on;
   t=[1:0.1:max(time)].';plot(t,xout(2)*exp(-t./xout(1))-xout(3),'k');
    legend('Measurement','Fit','Location','Best')
   title('T_2 Measurement')
   xlabel('TE [ms]')
   ylabel('Signal [a.u.]')
   grid on
   disp(['fitted T2: ' num2str(xout(1)) ' ms']);
   disp(['fitted M0: ' num2str(xout(2)*1e5) ' a.u.']);
   fitVal = xout(1);
elseif strcmp(mode,'FA')
    [est_A, ind] = max(specMax);
    
    if numel(voltage) >= ind*2
        est_omega = pi/voltage(ind*2);
    else
    	est_omega=pi/max(voltage);
    end
    
    % If signal is measured up to high voltages, assumption of sin-shape
    % doesn't suit anymore because of the B1-Inhomogeneity (see Collins2001
    % for example images)
    [m i] = findpeaks(specMax);
    if size(m,2) > 1 && m(2) < m(1)*0.9
        [M,Id] = min(specMax(i(1):i(2)));
        specMax = specMax(1:Id+i(1));
        voltage = voltage(1:Id+i(1));
    end

    xout=fminsearch('sin_abs',[est_omega est_A],[],voltage,specMax);

    figure,plot(voltage,specMax,'k+');hold on;
    t=[1:0.1:max(voltage)].';plot(t,abs(xout(2)*sin(xout(1).*t)),'k');
    legend('Measurement','Fit','Location','Best')

    title('FA Calibration')
    xlabel('Transmitter Voltage [V]')
    ylabel('Signal [a.u.]')
    grid on
    disp(['Reference Voltage: ' num2str(round(pi/xout(1))) ' V']);
    fitVal = round(pi/xout(1));
else
    disp('fit function not supported');
end
                        
        