% Generate lay file for 16-channel EEG from music project at HMRI MRI core 
% and convert from bipolar to monopolar montage referenced to O1 electrode
% (O1 electrode removed as a result)
global lay
cfg               = [];
cfg.elecfile      = 'MusicFavorite_Edit Markers_rem_sync_SSs_R_locs.sfp';
lay               = ft_prepare_layout(cfg);

% this is for matching position labels to data labels
newlabel = cell(length(lay.label),1);
for cp = 1:length(lay.label)
    newlabel{cp}  = regexprep(lay.label{cp},'-\w{2,3}','');
end
lay.label         = newlabel;
lay.cfg.channel   = newlabel;
% need to remove the O1 channel
o1IDX            = ismember(lay.label,'O1');
lay.cfg.channel   = lay.cfg.channel(~o1IDX);
lay.label         = lay.label(~o1IDX);
lay.pos           = lay.pos(~o1IDX,:);
lay.width         = lay.width(~o1IDX,:);
lay.height        = lay.height(~o1IDX,:);