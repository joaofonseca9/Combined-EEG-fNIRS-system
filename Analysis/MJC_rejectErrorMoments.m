function [EEG_out] = MJC_rejectErrorMoments(EEG_in, sub)
% For each subject, check if there is the need to remove any particular
% moment of the recording and remove it if so.

event_samp  = [EEG_in.event.latency];

% Sub-02.
if sub == "02"
    
    % Reject all up until the start of the automaticity test.
    % Get the moment to start cutting.
    start_cut1 = 1;
    % Get the moment to stop cutting.
    startTask_autodual_cued = event_samp(find(strcmp({EEG_in.event.type}, 's1706')==1));
    startTask_autodual_uncued = event_samp(find(strcmp({EEG_in.event.type}, 's1707')==1));
    end_cut1 = min(min(startTask_autodual_cued, startTask_autodual_uncued));
    
    % Identify when to remove moment at the end of the automaticity test up
    % until the beggining of the next task - in this case auto single.
    % Get the moment to start cutting.
    endTask_autodual_cued = event_samp(find(strcmp({EEG_in.event.type}, 's1714')==1));
    endTask_autodual_uncued = event_samp(find(strcmp({EEG_in.event.type}, 's1715')==1));
    start_cut2 = max(max(endTask_autodual_cued, endTask_autodual_uncued));
    % Get the moment to end cutting.
    startTask_autosingle_cued = event_samp(find(strcmp({EEG_in.event.type}, 's1702')==1));
    startTask_autosingle_uncued = event_samp(find(strcmp({EEG_in.event.type}, 's1703')==1));
    end_cut2 = min(min(startTask_autosingle_cued, startTask_autosingle_uncued));
    
    % Reject time between end of auto single and start of non-auto single.
    % Get the moment to start cutting.
    endTask_autosingle_cued = event_samp(find(strcmp({EEG_in.event.type}, 's1710')==1));
    endTask_autosingle_uncued = event_samp(find(strcmp({EEG_in.event.type}, 's1711')==1));
    start_cut3 = max(max(endTask_autosingle_cued, endTask_autosingle_uncued));
    % Get the moment to stop cutting.
    startTask_nonautosingle_cued = event_samp(find(strcmp({EEG_in.event.type}, 's1704')==1));
    startTask_nonautosingle_uncued = event_samp(find(strcmp({EEG_in.event.type}, 's1705')==1));
    end_cut3 = min(min(startTask_nonautosingle_cued, startTask_nonautosingle_uncued));
    
    % Identify when to remove moment at the start of the non-auto dual
    % task - eliminate first two trials out of 22.
    % Get the moment to start cutting.
    endTask_nonautosingle_cued = event_samp(find(strcmp({EEG_in.event.type}, 's1712')==1));
    endTask_nonautosingle_uncued = event_samp(find(strcmp({EEG_in.event.type}, 's1713')==1));
    start_cut4 = max(max(endTask_nonautosingle_cued, endTask_nonautosingle_uncued));
    % Get the moment to stop cutting.
    startTask_nonautodual_cued = event_samp(find(strcmp({EEG_in.event.type}, 's1708')==1));
    startTask_nonautodual_uncued = event_samp(find(strcmp({EEG_in.event.type}, 's1709')==1));
    startTask_nonautodual = sort([startTask_nonautodual_cued, startTask_nonautodual_uncued]);
    end_cut4 = startTask_nonautodual(3);
    
    % Reject time between end of non-auto dual task up until checkerboard
    % task.
    % Get the moment to start cutting.
    endTask_nonautodual_cued = event_samp(find(strcmp({EEG_in.event.type}, 's1716')==1));
    endTask_nonautodual_uncued = event_samp(find(strcmp({EEG_in.event.type}, 's1717')==1));
    start_cut5 = max(max(endTask_nonautodual_cued, endTask_nonautodual_uncued));
    % Get the moment to stop cutting.
    startTask_checkerboard = event_samp(find(strcmp({EEG_in.event.type}, 's1555')==1));
    end_cut5 = min(startTask_checkerboard);
    
    % Reject time from the end of the checkerboard up until the end of the
    % recording.
    % Get the moment to start cutting.
    endTask_checkerboard = event_samp(find(strcmp({EEG_in.event.type}, 's1500')==1));
    start_cut6 = min(endTask_checkerboard);
    % Get the moment to stop cutting.
    end_cut6 = size(EEG_in.data, 2);
    
    % Eliminate the moments from the signal.
    [EEG_aux] = eeg_eegrej(EEG_in,...
        [start_cut1 end_cut1-23*EEG_in.srate;...
        start_cut2+23*EEG_in.srate end_cut2-23*EEG_in.srate;...
        start_cut3+23*EEG_in.srate end_cut3-23*EEG_in.srate;...
        start_cut4+23*EEG_in.srate end_cut4-23*EEG_in.srate;...
        start_cut5+23*EEG_in.srate end_cut5-3*EEG_in.srate;...
        start_cut6+3*EEG_in.srate end_cut6]);  
    
    % Reject other bad moments of the data.
    EEG_aux = eeg_eegrej(EEG_aux, [49295 50778]);
    EEG_aux = eeg_eegrej(EEG_aux, [63349 64362]);
    EEG_aux = eeg_eegrej(EEG_aux, [81037 81993]);
    EEG_aux = eeg_eegrej(EEG_aux, [168513 169589]);
    EEG_aux = eeg_eegrej(EEG_aux, [169637 169996]);
    EEG_aux = eeg_eegrej(EEG_aux, [173668 174727]);
    EEG_aux = eeg_eegrej(EEG_aux, [187913 191166]);
    EEG_aux = eeg_eegrej(EEG_aux, [215832 218575]);
    EEG_aux = eeg_eegrej(EEG_aux, [259879 264275]);
    EEG_aux = eeg_eegrej(EEG_aux, [274310 275198]);
    EEG_aux = eeg_eegrej(EEG_aux, [300319 301366]);
    EEG_aux = eeg_eegrej(EEG_aux, [323766 328966]);
    EEG_aux = eeg_eegrej(EEG_aux, [328839 330007]);
    EEG_aux = eeg_eegrej(EEG_aux, [336533 338224]);
    EEG_aux = eeg_eegrej(EEG_aux, [369554 372199]);
    EEG_aux = eeg_eegrej(EEG_aux, [413710 416361]);
    EEG_aux = eeg_eegrej(EEG_aux, [457226 459342]);
    EEG_aux = eeg_eegrej(EEG_aux, [502212 503823]);
    EEG_aux = eeg_eegrej(EEG_aux, [541664 544936]);
    EEG_aux = eeg_eegrej(EEG_aux, [571840 574577]);
    EEG_aux = eeg_eegrej(EEG_aux, [577501 580256]);
    EEG_aux = eeg_eegrej(EEG_aux, [579620 580289; 580647 581979;...
        586563 588363; 590240 591110; 606444 607491; 609027 609581;...
        610179 611523; 618772 619271; 619823 621015; 632294 633790;...
        650191 652410; 656247 657317; 663512 664029; 674572 677175;...
        694563 695500; 703840 704625; 722803 724354; 737780 741758;...
        748543 750605; 758502 763196; 765507 769617; 788823 792934;...
        804686 806182; 822632 824456; 842693 843387; 859374 861600]);
    EEG_aux = eeg_eegrej(EEG_aux, [891687 892447; 973832 975401;...
        1017788 1019436; 1049360 1050065; 1050156 1051068;...
        1171115 1174076; 1196672 1198399; 1292194 1293665;...
        1302290 1303318; 1303379 1304990; 1386808 1388116;...
        1453204 1455867; 1515823 1516499; 1547972 1548666;...
        1548763 1549140; 1664044 1669484; 1675752 1676057;...
        1677315 1679498; 1687306 1688364; 1688832 1689440;...
        1707311 1707622; 1708485 1710084; 1712018 1713508;...
        1717685 1720044; 1720847 1721924; 1778079 1779205;...
        1819775 1820171; 1821058 1821782]);
    EEG_aux = eeg_eegrej(EEG_aux, [1808118 1808526; 1862340 1862869;...
        1863095 1865400; 1875068 1875859; 1907417 1907819;...
        1908586 1911584; 1926548 1927734; 1931336 1935190;...
        1956484 1957408; 1972418 1973256; 2013870 2017730;...
        2020157 2020734; 2021173 2022923; 2055875 2056332;...
        2085494 2089082; 2104331 2109201; 2110695 2111801;...
        2133976 2135126; 2187780 2188778; 2201979 2204417;...
        2236616 2241029; 2241565 2242781; 2274115 2277192;...
        2277339 2279680; 2298727 2301610; 2316760 2318324;...
        2332662 2334388; 2347037 2348934; 2363718 2366771;...
        2405603 2407330; 2410230 2411593; 2434940 2435615;...
        2452236 2458278; 2472710 2473519; 2510873 2511742;...
        2512052 2514242; 2523454 2525862; 2724891 2725731;...
        2751306 2752303; 2798111 2799145; 2844610 2845499;...
        2848508 2850594; 2892382 2894493; 2972241 2972673;...
        2973170 2974398; 3148326 3148655; 3149118 3151021;...
        3190704 3191556; 3329200 3329583; 3377391 3379538;...
        3487296 3489868]);
    EEG_out = EEG_aux;
    
% Sub-76.
elseif sub == "76"
    
    % Reject all up until the start of the second trial of the automaticity
    % test.
    % Get the moment to start cutting.
    start_cut1 = 1;
    % Get the moment to stop cutting.
    startTask_autodual_cued = event_samp(find(strcmp({EEG_in.event.type}, 's1706')==1));
    startTask_autodual_uncued = event_samp(find(strcmp({EEG_in.event.type}, 's1707')==1));
    startTask_autodual = sort([startTask_autodual_cued, startTask_autodual_uncued]);
    end_cut1 = startTask_autodual(2);
    
    % Identify when to remove moment at the end of the automaticity test up
    % until the beggining of the next task - in this case non-auto single.
    % Get the moment to start cutting.
    endTask_autodual_cued = event_samp(find(strcmp({EEG_in.event.type}, 's1714')==1));
    endTask_autodual_uncued = event_samp(find(strcmp({EEG_in.event.type}, 's1715')==1));
    start_cut2 = max(max(endTask_autodual_cued, endTask_autodual_uncued));
    % Get the moment to end cutting.
    startTask_autosingle_cued = event_samp(find(strcmp({EEG_in.event.type}, 's1704')==1));
    startTask_autosingle_uncued = event_samp(find(strcmp({EEG_in.event.type}, 's1705')==1));
    end_cut2 = min(min(startTask_autosingle_cued, startTask_autosingle_uncued));
    
    % Reject time between end of non-auto single and start of non-auto dual.
    % Get the moment to start cutting.
    endTask_nonautosingle_cued = event_samp(find(strcmp({EEG_in.event.type}, 's1712')==1));
    endTask_nonautosingle_uncued = event_samp(find(strcmp({EEG_in.event.type}, 's1713')==1));
    start_cut3 = max(max(endTask_nonautosingle_cued, endTask_nonautosingle_uncued));
    % Get the moment to stop cutting.
    startTask_nonautodual_cued = event_samp(find(strcmp({EEG_in.event.type}, 's1708')==1));
    startTask_nonautodual_uncued = event_samp(find(strcmp({EEG_in.event.type}, 's1709')==1));
    end_cut3 = min(min(startTask_nonautodual_cued, startTask_nonautodual_uncued));
    
    % Reject time between end of non-auto dual task up until the start of
    % the auto single.
    % Get the moment to start cutting.
    endTask_nonautodual_cued = event_samp(find(strcmp({EEG_in.event.type}, 's1716')==1));
    endTask_nonautodual_uncued = event_samp(find(strcmp({EEG_in.event.type}, 's1717')==1));
    start_cut4 = max(max(endTask_nonautodual_cued, endTask_nonautodual_uncued));
    % Get the moment to end cutting.
    startTask_autosingle_cued = event_samp(find(strcmp({EEG_in.event.type}, 's1702')==1));
    startTask_autosingle_uncued = event_samp(find(strcmp({EEG_in.event.type}, 's1703')==1));
    end_cut4 = min(min(startTask_autosingle_cued, startTask_autosingle_uncued));
    
    % Reject time between end of auto single and start of checkerboard.
    % Get the moment to start cutting.
    endTask_autosingle_cued = event_samp(find(strcmp({EEG_in.event.type}, 's1710')==1));
    endTask_autosingle_uncued = event_samp(find(strcmp({EEG_in.event.type}, 's1711')==1));
    start_cut5 = max(max(endTask_autosingle_cued, endTask_autosingle_uncued));
    % Get the moment to stop cutting.
    startTask_checkerboard = event_samp(find(strcmp({EEG_in.event.type}, 's1555')==1));
    end_cut5 = min(startTask_checkerboard);
    
    % Reject time between end of checkerboard and end of recording.
    % Get the moment to start cutting.
    endTask_checkerboard = event_samp(find(strcmp({EEG_in.event.type}, 's1500')==1));
    start_cut6 = min(endTask_checkerboard);
    % Get the moment to stop cutting.
    end_cut6 = size(EEG_in.data, 2);

[EEG_out] = eeg_eegrej(EEG_in,...
    [start_cut1 end_cut1-23*EEG_in.srate;...
    start_cut2+23*EEG_in.srate end_cut2-23*EEG_in.srate;...
    start_cut3+23*EEG_in.srate end_cut3-23*EEG_in.srate;...
    start_cut4+23*EEG_in.srate end_cut4-23*EEG_in.srate;...
    start_cut5+23*EEG_in.srate end_cut5-3*EEG_in.srate;...
    start_cut6+3*EEG_in.srate end_cut6]);

else
    
    EEG_out = EEG_in;

end

StartBlock_AutomaticSequence_Cued = length(find(strcmp({EEG_out.event.type}, 's1702')==1))
StartBlock_AutomaticSequence_Uncued = length(find(strcmp({EEG_out.event.type}, 's1703')==1))
StartBlock_NonAutomaticSequence_Cued = length(find(strcmp({EEG_out.event.type}, 's1704')==1))
StartBlock_NonAutomaticSequence_Uncued = length(find(strcmp({EEG_out.event.type}, 's1705')==1))

StartBlock_AutomaticSequence_Dual_Cued = length(find(strcmp({EEG_out.event.type}, 's1706')==1))
StartBlock_AutomaticSequence_Dual_Uncued = length(find(strcmp({EEG_out.event.type}, 's1707')==1))
StartBlock_NonAutomaticSequence_Dual_Cued = length(find(strcmp({EEG_out.event.type}, 's1708')==1))
StartBlock_NonAutomaticSequence_Dual_Uncued = length(find(strcmp({EEG_out.event.type}, 's1709')==1))

EndBlock_AutomaticSequence_Cued = length(find(strcmp({EEG_out.event.type}, 's1710')==1))
EndBlock_AutomaticSequence_Uncued = length(find(strcmp({EEG_out.event.type}, 's1711')==1))
EndBlock_NonAutomaticSequence_Cued = length(find(strcmp({EEG_out.event.type}, 's1712')==1))
EndBlock_NonAutomaticSequence_Uncued = length(find(strcmp({EEG_out.event.type}, 's1713')==1))

EndBlock_AutomaticSequence_Dual_Cued = length(find(strcmp({EEG_out.event.type}, 's1714')==1))
EndBlock_AutomaticSequence_Dual_Uncued = length(find(strcmp({EEG_out.event.type}, 's1715')==1))
EndBlock_NonAutomaticSequence_Dual_Cued = length(find(strcmp({EEG_out.event.type}, 's1716')==1))
EndBlock_NonAutomaticSequence_Dual_Uncued = length(find(strcmp({EEG_out.event.type}, 's1717')==1))

StartBlock_Checkerboard = length(find(strcmp({EEG_out.event.type}, 's1555')==1))
EndBlock_Checkerboard = length(find(strcmp({EEG_out.event.type}, 's1500')==1))

end