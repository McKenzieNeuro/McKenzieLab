function reSave_sparcNet(fil)

dirN = fileparts(fil);
Fs=200;ww=2*Fs;

    disp(['add #',num2str(i),' ',files{i,1}])
    
    mat_file=matfile([data_dir,files{i,1}]);
    [xx,num_points]=size(mat_file,'data');
    num_seg=ceil(num_points/ww);

    featue_file=[sn_dir,strrep(files{i,1},'.mat','_score.csv')]; 
    Y=cell2mat(table2cell(readtable(featue_file)));
    Y=[repmat(Y(1,:),2,1);Y];
    Y=[Y;repmat(Y(end,:),num_seg-size(Y,1),1)];
    save([res_dir,files{i,1}],'Y')

