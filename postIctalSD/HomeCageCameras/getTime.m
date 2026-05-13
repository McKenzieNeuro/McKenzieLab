function timeVec = getTime()
    Now = clock;
    timeVec = (Now(1,4)*3600)+(Now(1,5)*60)+(Now(1,6));
end
