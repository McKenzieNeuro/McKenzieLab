function checkMatprefs()

    v = version;
    v = symsepchar(v,' ');
    v = v{1,2}(1,2:end-1);

    user = compiler.UserInfo();
    user = user.UserID;

    path = ['C:\Users\',user,'\AppData\Roaming\MathWorks\MATLAB\',v];
    
    try
        load([path,'\matlabprefs.mat'])
    catch
        eval(['copyfile ',path,'\matlabprefs_2.mat ',path,'\matlabprefs.mat'])
    end
end