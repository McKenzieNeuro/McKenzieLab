function data = sm_LoadPickle(fname)
% in your python environment pip install ruamel.yaml
% and before you need to set your python env in matlab 
% % e.g. pyversion('C:\Users\samckenzie\Anaconda3\python.exe')


fid = py.open(fname,'rb');
data = py.pickle.load(fid);
end