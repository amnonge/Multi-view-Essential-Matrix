
clear

%Uncomment this to run on all data sets

% dataList={  'Vienna_Cathedral', 'Piazza_del_Popolo','NYC_Library' ,'Alamo','Madrid_Metropolis' , 'Yorkminster',...
%     'Montreal_Notre_Dame','Tower_of_London','Ellis_Island','Notre_Dame'};

dataList={  'Ellis_Island'};

% Define the parameters of the model
ParamOpt = struct; 
ParamOpt.col = 0.985; %collinear threshold, set to 1 to include all cameras in the scene
ParamOpt.rot = 1.1; %rotation consistency threshold, set to 3 to include all cameras
ParamOpt.tran =  1;


DataRange = [1:length(dataList)];
Results_for_table = struct;
%write here the path to your data
path2 = '1DSFMEdata';

for iii = DataRange
    
    dataset = dataList{iii};
    %run the pipline on the data set
    results=runPipline(dataset,ParamOpt,path2);
    
    %save the results into a struct to be used later.
    Results_for_table(iii).name = dataset;
    Results_for_table(iii).optTime=results.optTime;
    Results_for_table(iii).finalTime=results.finalTime;
    Results_for_table(iii).Rotation_frob_mean = results.RF1;
    Results_for_table(iii).Rotation_frob_median = results.RF2;
    Results_for_table(iii).Translation_mean = results.T1;
    Results_for_table(iii).Translation_median = results.T2;
    Results_for_table(iii).Rotation_deg_mean = results.RD1;
    Results_for_table(iii).Rotation_deg_median = results.RD2;
end

t = table({Results_for_table.name}',[Results_for_table.Rotation_frob_mean]',[Results_for_table.Rotation_deg_mean]',[Results_for_table.Translation_mean]',[Results_for_table.Translation_median]'...
    ,[Results_for_table.optTime]',[Results_for_table.finalTime]','VariableNames',{'Dataset','Rotation_frob_mean','Rotation_deg_mean','Translation_mean','Translation_median','optTime','finalTime'});
writetable(t,'table.csv');
