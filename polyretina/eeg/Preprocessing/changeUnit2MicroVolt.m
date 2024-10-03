function EEG = changeUnit2MicroVolt(EEG, cfg)
switch lower(cfg.recording_unit)
    case 'volt'
        coeff = 10^6;
    otherwise
        error('Unit not recognized')
end

EEG.data = EEG.data.*coeff;
if ~isempty(EEG.icaact)
    EEG.icaact = EEG.icaact.*coeff;
end
end