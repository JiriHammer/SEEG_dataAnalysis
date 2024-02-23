function mytitle = textOnFig_chnlInfo(params, plotInfo)
% returns text (string) of channel info
% placed on top of figure

% (c) Jiri, Jun20

if ~isfield(plotInfo, 'channel')
    H_ch = struct;
else
    H_ch = plotInfo.channel;
end

%% defaults
if ~isfield(H_ch, 'neurologyLabel') 
    H_ch.neurologyLabel = 'n.a.'; 
end

%% text string (first row): channel name, MNIs, anatomy assignemnts
mytitle = [params.storage.subjTag ', ch = ' H_ch.name '(' num2str(plotInfo.thisCh) ')'];

% MNI
if isfield(H_ch, 'MNI_x')
    if ~isempty(H_ch.MNI_x)
        mni = [H_ch.MNI_x, H_ch.MNI_y, H_ch.MNI_z];
        mytitle = [mytitle ', MNI = [' num2str(mni,'% 3.0f') ']'];
    end
end

% anatomic assignment
if isfield(H_ch, 'ass_brainAtlas')
    anatStr = shortenChannelName(H_ch.ass_brainAtlas, 6);
%     mytitle = [mytitle ', anat = ' H_ch.ass_brainAtlas];
    mytitle = [mytitle ', anat = ' anatStr];
end

% cytoarchitectonic map
if isfield(H_ch, 'ass_cytoarchMap')
    mytitle = [mytitle ', cyto = ' H_ch.ass_cytoarchMap];
end

% neurologist
if isfield(H_ch, 'neurologyLabel')
    mytitle = [mytitle ', nrlg = ' H_ch.neurologyLabel];
end

%% anatomic assignments from atlas (second row)
mytitle_ass = [];

% Yeo7
if isfield(H_ch, 'ass_yeo7_name')
%     mytitle = [mytitle ', yeo7 = ' H_ch.ass_yeo7_name ' (d=' num2str(H_ch.ass_yeo7_dist,2) ' mm)'];
    mytitle_ass = [mytitle_ass 'yeo7 = ' H_ch.ass_yeo7_name ' (d=' num2str(H_ch.ass_yeo7_dist,2) ' mm)'];
end

% Yeo17
if isfield(H_ch, 'ass_yeo17_name')
%     mytitle = [mytitle ', yeo17 = ' H_ch.ass_yeo17_name ' (d=' num2str(H_ch.ass_yeo17_dist,2) ' mm)'];
    mytitle_ass = [mytitle_ass ', yeo17 = ' H_ch.ass_yeo17_name ' (d=' num2str(H_ch.ass_yeo17_dist,2) ' mm)'];
end

% Mars
if isfield(H_ch, 'ass_mars_name')
%     mytitle = [mytitle ', yeo7 = ' H_ch.ass_mars_name ' (d=' num2str(H_ch.ass_mars_dist,2) ' mm)'];
    mytitle_ass = [mytitle_ass ', mars = ' H_ch.ass_mars_name ' (d=' num2str(H_ch.ass_mars_dist,2) ' mm)'];
end

% GM / WM / CSF / BONE
if isfield(H_ch, 'ass_gmwm_name')
    gw_str = [];
    for gw = 1:size(H_ch.ass_gmwm_name,2)
        gw_str = [gw_str ', ' H_ch.ass_gmwm_name{gw} '(' num2str(round(H_ch.ass_gmwm_dist(gw)*100)) '%)'];
    end
    mytitle_ass = [mytitle_ass gw_str];
end

%% print text
if nargout == 0
    tx = axes('visible','off', 'position',[0 0 1 1]);
    mytitle = strrep(mytitle, '_','\_');
    text(plotInfo.txt_pos(1), plotInfo.txt_pos(2), mytitle, 'fontsize', 16, 'fontw', 'bold'); 
    
    % second row of text
    if ~isempty(mytitle_ass)
        tx = axes('visible','off', 'position',[0 0 1 1]);
        mytitle = strrep(mytitle_ass, '_','\_');
        text(plotInfo.txt_pos(1), plotInfo.txt_pos(2)-0.03, mytitle, 'fontsize', 16, 'fontw', 'bold'); 
    end 
end
