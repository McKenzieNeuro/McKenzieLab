function setStim(TCP,enable,config)

    % Set up a channel's stim parameters
    % TCP = Intan TCP object
    % enable = turn stim on or off, 0 = off, 1 = on
    % config = structure with fields:   
        % port = port on Intan: 'a','b','c','d'
        % fKey = F key to trigger stim
        % chan = channel number (base zero), [c1,c2] for bipolar
        % amp = Stim amplitude in uA
        % polarity = leading current, -1 = cathode first, 1 = anode first
        % duration = pulse duration in uS
        % pulseNum = number of pulses
        % period = pulse train period in mS
    
    for i = 1:size(config,2)
        for ii = 1:length(config(i).chan)
            if ~isnan(config(i).chan(ii))

                if config(i).chan(ii) < 10
                    chanStr = [config(i).port,'-00',num2str(config(i).chan(ii))];
                else 
                    chanStr = [config(i).port,'-0',num2str(config(i).chan(ii))];
                end

                if enable == 1              

                    write(TCP,uint8(['set ',chanStr,'.stimenabled true;']));
                    write(TCP,uint8(['set ',chanStr,'.source KeyPressF',num2str(config(i).fKey),';']));
                    write(TCP,uint8(['set ',chanStr,'.FirstPhaseAmplitudeMicroAmps ',num2str(config(i).amp),';']));
                    write(TCP,uint8(['set ',chanStr,'.SecondPhaseAmplitudeMicroAmps ',num2str(config(i).amp),';']));
    
                    if config(i).polarity(ii) == 1
                        write(TCP,uint8(['set ',chanStr,'.Polarity PositiveFirst;']));
                    elseif config(i).polarity(ii) == -1
                        write(TCP,uint8(['set ',chanStr,'.Polarity NegativeFirst;']));
                    end
    
                    write(TCP,uint8(['set ',chanStr,'.FirstPhaseDurationMicroseconds ',num2str(round(config(i).duration/2)),';']));
                    write(TCP,uint8(['set ',chanStr,'.SecondPhaseDurationMicroseconds ',num2str(round(config(i).duration/2)),';']));

                    if config(i).pulseNum == 1
                        write(TCP,uint8(['set ',chanStr,'.PulseOrTrain SinglePulse;']));
                    else
                        write(TCP,uint8(['set ',chanStr,'.PulseOrTrain PulseTrain;']));
                        write(TCP,uint8(['set ',chanStr,'.NumberOfStimPulses ',num2str(config(i).pulseNum),';']));
                        write(TCP,uint8(['set ',chanStr,'.PulseTrainPeriodMicroseconds ',num2str(config(i).period*1000),';']));
                    end
                elseif enable == 0
                    write(TCP,uint8(['set ',chanStr,'.stimenabled false;']));
                end
            end
        end
    end
end