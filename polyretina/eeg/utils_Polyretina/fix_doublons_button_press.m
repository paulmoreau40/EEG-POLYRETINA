function [EEG_cleaned] = fix_doublons_button_press(EEG)
% Custom function to remove button press doublons the data according to your needs
%
% Inputs:
% EEG           - EEG struct where Button Press doublons need to be removed
%
% Outputs:
% EEG_cleaned   - EEG struct where no more doublons

% Getting indices where button press tag 
button_press_mask = strcmp({EEG.event.type}, 'Button_Press');

doublon_count = 0;

for j = 1:length(button_press_mask)-1
    
    % Looping to see if 2 consecutive button presses
    if button_press_mask(j) % We are at a button press
        if button_press_mask(j+1) % The next index is also a button press: we are in a doublon scenario
            button_press_mask(j+1) = 0;
            doublon_count = doublon_count + 1;
        else % We only have a single button press
            button_press_mask(j) = 0;
        end
    end
    
end

% Removing Doublons
disp(['Found ' num2str(doublon_count) ' Button Press doublons, removing them...'])

EEG_cleaned = pop_editeventvals(EEG, 'delete', find(button_press_mask));


end

