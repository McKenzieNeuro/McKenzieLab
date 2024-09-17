function tim = sm_timeFromEvent(ts,ev)

tim = nan(size(ts));

for i = 1:length(ev)-1
    kp = ts>=ev(i) & ts<ev(i+1);
    tim(kp) = ts(kp) - ev(i);
end

kp = ts>=ev(end);
tim(kp) = ts(kp) - ev(end);


end