clear
datasets={  'Vienna_Cathedral', 'Piazza_del_Popolo','NYC_Library' ,'Alamo','Madrid_Metropolis' , 'Yorkminster',...
    'Montreal_Notre_Dame','Tower_of_London','Ellis_Island','Notre_Dame'};


means=[];
medians=[];

meansBA=[];
mediansBA=[];
Nrs=[];
for i=1:length(datasets)
 [RF1,RF2,T1,T2,RD1,RD2,N]=GetFinalResults(datasets{i});
 means=[means;T1];
 medians=[medians;T2];
 Nrs=[Nrs;N];
end

T=table(datasets',means,medians,Nrs);
writetable(T,'tableFinal.csv');
