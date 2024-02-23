function [EC,F,T]=MVAR_v2(signal, time, labels, setting)
% INPUTS===================================================================
% signal... EEG (time X channel X trials)
% time... time vector in seconds (time X 1 X 1)
% labels... (1 x trials) defines trials to avaraging
% setting... structure of individual settings
% setting.ww... segmentation window (sec.)
% setting.nn... negative overlap x100% 0.9~90%
% setting.bandwidth... [bf1_min bf1_max; bf2_min bf2_max;...] bands of
%                        partial MVAR estimation
% setting.pooled... 0 or 1; parallel computing (parfor) no/yes 
% setting.ec_matrices={'SdDTF','DTF','ffDTF','PDC','COH'}; select metrices
% to output EC
% OUTPUTS: 
% EC... effective connectivity structure of metrices (used abs value, can be modified in code)
%             (CH x CH x F x T x label) 
% METRICES by the "mvfreqz_faster.m", see details in function
%   S   	power spectrum
%   h	transfer functions, abs(h.^2) is the non-normalized DTF [11]
%   PDC 	partial directed coherence [2]
%   DC  	directed coupling	
%   COH 	coherency (complex coherence) [5]
%   DTF 	directed transfer function
%   pCOH 	partial coherence
%   dDTF 	direct Directed Transfer function (SdDTF)
%   ffDTF full frequency Directed Transfer Function 
%   pCOH2 partial coherence - alternative method 
%   GGC	a modified version of Geweke's Granger Causality [Geweke 1982]
%	   !!! it uses a Multivariate AR model, and computes the bivariate GGC as in [Bressler et al 2007]. 
%	   This is not the same as using bivariate AR models and GGC as in [Bressler et al 2007]
%   Af	Frequency transform of A(z), abs(Af.^2) is the non-normalized PDC [11]
%   PDCF 	Partial Directed Coherence Factor [2]
%   GPDC 	Generalized Partial Directed Coherence [9,10]
%
% F... frequency vector of computed freq.
% T... time-segments indexes at the start of window

if ~isfield(setting, 'doFreqBandFiltering')
    if size(setting.bandwidth,1) > 1
        setting.doFreqBandFiltering = true;    % added by JH (Nov22)
    else
        setting.doFreqBandFiltering = false;
    end
end
warning('off','MATLAB:nearlySingularMatrix');

signal=single(signal); % in single precesion to memory save
fs=1/median(diff(time)); % sampling frequency from time vector
disp(['fs=' num2str(fs) 'Hz'])
ww=round(setting.ww*fs);  % window size (samples)
nn=round(setting.nn*ww); % noverlap (samples)

index=1:ww-nn:size(signal,1)-ww+1;
t_index=time(index);
bandwidth=setting.bandwidth;

T=time(index); % time of segments (start of window)
T=T(:);


fs_new=fs;
unique_labels=unique(labels);

% memory reservation of effective connectivity (EC)
F=[];
for i=1:size(bandwidth,1)
    F=[F bandwidth(i,1):bandwidth(i,2)]; % generate frequency axis, that defines also size of EC matrices
end

% reservation of memory by the required EC mtrices
for i=1:length(setting.ec_matrices)
    % dynamic structure reference   
    EC.(setting.ec_matrices{i})=single(NaN(size(signal,2),size(signal,2),length(F),length(index),length(unique_labels)));
end
% possible outputs of mvfreqz_faster.m function in order
ec_out={'S','h','PDC','COH','DTF','DC','pCOH','dDTF','ffDTF', 'pCOH2', 'PDCF', 'coh','GGC','Af','GPDC'};
for i=1:length(setting.ec_matrices)
    if sum(strcmp(setting.ec_matrices{i},ec_out))==0
        error(['Required EC matrices "' setting.ec_matrices{i} '" is not in supported list. See HELP MVAR_v2'])
    end
end
[~,used_ec]=intersect(ec_out,setting.ec_matrices);

disp('=== MVAR ===')
for bw=size(bandwidth,1):-1:1 % MVAR estimation for narrow-bands from higher to lower
    % the MVAR estimation is divided to sub-bands defined by setting.bandwidt. 
    % Algorithm starts in upper bandwidth that band-pass filtering is possible 
    % without signal decimation. In next itteration of lower nabds, if the bandwidth ends
    % above (fs/2)/3, signal is decimated to fs_new and band-pass filter is
    % redesigned. Change of fs shift indexes of time-segmentation,
    % therefore must be found new indexes for new timeline "time_dec".
    % Time segments of all trials are serialized in individual channels.
    % MVAR estimates transfer matrices by trials in same labels, which replace
    % averaging and save huge computational time.  
    
    % decimate the signal for band-pass filter design and separate MVAR
    % computation
    
    disp(['sub-band: ' num2str(bandwidth(bw,1)) '-' num2str(bandwidth(bw,2)) ' Hz'])
    
    if setting.doFreqBandFiltering
        if 3*bandwidth(bw,2)<fs_new/2
            r=fs_new/(3*bandwidth(bw,2));
            dec_fs=round(fs_new/r);
        else
            dec_fs=fs_new;
        end

        if fs_new~=dec_fs
            signal_dec=[];
            for TR=1:size(signal,3)
                signal_dec(:,:,TR)=single(high_order_resample(double(signal(:,:,TR)),fs_new,dec_fs));
            end
            signal=signal_dec;
            fs_new=dec_fs;
        end
        
        % filter designer
        Wp=2*bandwidth(bw,:)/fs_new;
        if 0.1*fs_new>bandwidth(bw,1)
            Wp=2*bandwidth(bw,2)/fs_new;
            Ws=(bandwidth(bw,2)+0.1*fs_new)/fs_new;
            warning('band-pass is so wide, lower frequency is ignored, use low-pass only ')
        else
            Ws=2*[bandwidth(bw,1)-0.1*fs_new,  bandwidth(bw,2)+0.1*fs_new]/fs_new; Ws(Ws<=0)=0.01; Ws(Ws>=1)=0.99;
        end
        Rp=3;
        Rs=60;
        [n,Wn] = buttord(Wp,Ws,Rp,Rs);
        [bp,ap] = butter(n,Wn);
        % filter test
        [h,~]=freqz(bp,ap,1000,fs_new);
        if sum(isnan(h))>0 || sum(abs(h)>(1+1e-6))
            error('wrong filter design')
        end

        % signal filtration in narrow-band
        signal_f=single(zeros(size(signal)));
        for TR=1:size(signal,3)
            signal_f(:,:,TR)=single(filtfilt(bp,ap,double(signal(:,:,TR))));
        end 
    else
        signal_f = signal;   
    end
    % ---------------------------------------
    
    
    % MVAR model order
    AR_r=diff(bandwidth(bw,:)); 
    if AR_r>setting.order_max 
        AR_r=setting.order_max;     % maximal order
    end  
    
    % window segmentation should not agree with time samples in original signal fs,
    % therefore indexes of time segment must be rounded and reindexed 
    time_dec=linspace(time(1),time(end),size(signal,1));
    
    % time segmentation: MVAR for each time-window. The new proces serialize time-window from same labels epochs. 
    % Averaging is processed during MVAR estimation (I hope). 
    
    f=bandwidth(bw,1):1:bandwidth(bw,end); % freqencies to compute effective connectivity range
    
    [~,fidx]=intersect(F,f);
    switch setting.pooled
        case 0 % without pooling
            % dynamic definition of outputs by the settings of
            % required connectivity matrices
            out_code='[';
            for i=1:max(used_ec)
                if isempty(intersect(i,used_ec))
                    out_code=[out_code '~,'];
                else
                    out_code=[out_code 'fvar.' ec_out{i} ','];
                end
            end
            out_code(end)=']'; % generate code-call of ouputs for eval
            
            
            for L=1:length(unique_labels) % separately for labels type
                
                SEG={}; % select all time-segments and test correlation
                for ti=1:length(index)
                    start=find(time_dec>=t_index(ti),1);
                    stop=find(time_dec<=(t_index(ti)+(setting.ww)-1/fs),1,'last');
                    seg=signal_f(start:stop,:,labels==unique_labels(L)); % (T x CH x trial)
                    
                    % detrending in each window (added by Jiri, Nov23)
                    for tr = 1:size(seg,3)
                        seg(:,:,tr) = detrend(seg(:,:,tr)); % detrend
                    end
%                     seg = zscore(seg,0,1);
                    
                    % equal same segments from trials type are serialized to SEG{ti} 
                    SEG{ti,1}=reshape(permute(seg,[2 1 3]),[size(seg,2) size(seg,1)*size(seg,3)])'; 
                    % [trial 1: T x CH; trial 2: T x CH; ...]
                    % trials serialization. Trials in 3rd dim. are putted bellow
                    % 1st
                    % example:
                    % A=repmat(1:10,[100 1]);
                    % A=cat(3,A,10*A,100*A);
                    % A=reshape(permute(A,[2 1 3]),[size(A,2) size(A,1)*size(A,3)])';
                    
                    % test of correlated channels and adding noise
                    CC=corr(SEG{ti,1});
                    CC(logical(eye(size(CC,1))))=0;
                    [r,c]=find(CC>0.99);
                    if ~isempty(r)
                        CR=sort([r,c],2); CR=unique(CR,'rows');
                        warning('channels are highly correlated, remove on of pair. The noise component were added.')
                        disp(CR)
                        % (un)commented by Jiri: Nov22
%                         SEG{ti,1}(:,CR(:,1))=SEG{ti,1}(:,CR(:,1))+0.05*repmat(std(SEG{ti,1}(:,CR(:,1)),[],1),[size(SEG{ti,1},1) 1])...
%                             .*randn(size(SEG{ti,1},1),size(CR,1)); % add 5% of noise
                    end
                end
                
                for ti=1:length(index) % for each time segment
                    start=find(time_dec>=t_index(ti),1);
                    stop=find(time_dec<=(t_index(ti)+(setting.ww)-1/fs),1,'last');
                                        
                    % === MVAR ===
                    [AR,~,PE]=mvar_simplified(SEG{ti,1},AR_r); % multi-chanel AR model
                    
                    % === effective connectivity ===
                    %  [~,~, fvar_PDC, fvar_COH, fvar_DTF, ~, ~, fvar_SdDTF, fvar_ffDTF] = mvfreqz_faster(...
                    %     eye(size(SEG{ti,1},2)),[eye(size(SEG{ti,1},2)),-AR],PE(:,AR_r*size(SEG{ti,1},2)+(1:size(SEG{ti,1},2))),f,fs_new); % CH x CH x f
                    eval([out_code '=mvfreqz_faster(eye(size(SEG{ti,1},2)),[eye(size(SEG{ti,1},2)),-AR],PE(:,AR_r*size(SEG{ti,1},2)+(1:size(SEG{ti,1},2))),f,fs_new);']);
                    % =====================================================
                    
                    
                    % store results
                    for i=1:length(setting.ec_matrices)
                        EC.(setting.ec_matrices{i})(:,:,fidx,ti,L)=single(abs(fvar.(setting.ec_matrices{i})));
                    end
%                     if isfield(EC,'DTF');   EC.DTF(:,:,fidx,ti,L)=single(fvar.DTF); end
%                     ... etc.
                end
            end
            
        case 1 % similar as "case 0" but rewritted for parpool 
            % dynamic definition of outputs by the settings of
            % required connectivity matrices
            out_code='[';
            for i=1:max(used_ec)
                if isempty(intersect(i,used_ec))
                    out_code=[out_code '~,'];
                else
                    out_code=[out_code 'fvar.' ec_out{i} '(:,:,:,ti),'];
                end
            end
            out_code(end)=']';
                
            for L=1:length(unique_labels)
                % reservation of memory by the required EC mtrices
                for i=1:length(setting.ec_matrices)
                    % dynamic structure reference
                    fvar.(setting.ec_matrices{i})=single(zeros(size(signal,2),size(signal,2),length(f),length(index)));
                end
%                 fvar.PDC=single(zeros(size(signal,2),size(signal,2),length(f),length(index)));
%                 ... etc.
                
                SEG={}; % select all time-segments sended to parloop
                for ti=1:length(index)
                    start=find(time_dec>=t_index(ti),1);
                    stop=find(time_dec<=(t_index(ti)+(setting.ww)-1/fs),1,'last');
                    seg=signal_f(start:stop,:,labels==unique_labels(L)); % (T x CH x trial)
                    % trials serialization. Trials in 3rd dim. are putted bellow
                    % 1st
                    % example:
                    % L=repmat(1:10,[100 1]);
                    % L=cat(3,L,10*L,100*L);
                    % L=reshape(permute(L,[2 1 3]),[size(L,2) size(L,1)*size(L,3)])';
                    SEG{ti,1}=reshape(permute(seg,[2 1 3]),[size(seg,2) size(seg,1)*size(seg,3)])'; % [trial 1: T x CH;
                                                                                                    %  trial 2: T x CH;...]                                                                                                    
                    % test of correlated channels and adding noise
                    CC=corr(SEG{ti,1});
                    CC(logical(eye(size(CC,1))))=0;
                    [r,c]=find(CC>0.99);
                    if ~isempty(r)
                        CR=sort([r,c],2); CR=unique(CR,'rows');
                        warning('channels are highly correlated, remove on of pair. The noise component were added.')
                        disp(CR)
                        
                        SEG{ti,1}(:,CR(:,1))=SEG{ti,1}(:,CR(:,1))+0.05*repmat(std(SEG{ti,1}(:,CR(:,1)),[],1),[size(SEG{ti,1},1) 1])...
                            .*randn(size(SEG{ti,1},1),size(CR,1)); % +5% šumu
                    end
                end
                
                % === MVAR === in parallel pooling (by time-segments ti)
                AR=zeros(size(SEG{1},2),size(SEG{1},2)*AR_r,length(index));
                PE=zeros(size(SEG{1},2),size(SEG{1},2)*(AR_r+1),length(index));
                parfor ti=1:length(index)
                    warning('off','MATLAB:nearlySingularMatrix');
                    
                    [AR(:,:,ti),~,PE(:,:,ti)]=mvar_simplified(SEG{ti},AR_r); % multi-chanel AR model
                    % === effective connectivity ===
                end
                % compute EC only by the required metrices
                for ti=1:length(index)
                    % dynamic definition of outputs ba the settings of
                    % required connectivity matrices
                    eval([out_code '=mvfreqz_faster(eye(size(SEG{ti,1},2)),[eye(size(SEG{ti,1},2)),-AR(:,:,ti)],PE(:,AR_r*size(SEG{ti,1},2)+(1:size(SEG{ti,1},2)),ti),f,fs_new);']);
                    
                    % =====================================================
                end
                % strore result
                for i=1:length(setting.ec_matrices)
                    EC.(setting.ec_matrices{i})(:,:,fidx,:,L)=single(abs(fvar.(setting.ec_matrices{i})));
                end
                % EC.DTF(:,:,fidx,:,L)=single(fvar.DTF);
                % ... etc.
            end
    end
end




