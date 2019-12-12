function results=runPipline(dataset,ParamOpt,path)
    load(sprintf([path '/%s.mat'], dataset));
    InitTime = tic;
    
    %Building the triplet cover and removing noisy edges
    [finalTriplets,v,~,tts] =  buildTriplets(pointMatchesInliers1,EN,Hmat,TijMat,ParamOpt);
    nodesNum=length(unique(finalTriplets));


    graphTime = toc(InitTime);
    initOpt = tic;
    % Main optimization function
    [Xs,Y,Xs1] = Optimization( EN,finalTriplets(1:end,1:3),num2cell(tts));
    
    optTime = toc(initOpt);
    % Extracting Rotation and translations from a set of consistent 3-view
    % matrices
    [ERR,ourR,ourT,camList,Fnew] = constructE(Xs,v,finalTriplets(:,1:3));

    finalTime = toc(InitTime);



    % Post processing, comparing to the ground truth
    [~,M,R_gt,T_gt,K_gt] = getSubGraph(M,EN,R_gt,T_gt,K_gt,camList);
    namesList =   namesList(camList);
    [RF1,RF2,T1,T2,RD1,RD2,scale]= getBestError(R_gt,ourR,T_gt,ourT ,1:length(R_gt));
    
    save(sprintf('results/%s',dataset),'R_gt','ourR','T_gt','ourT','namesList');
    
    % Display errors
    display('R frob mean|R frob median|T mean|T median| R degrees mean|R degrees median')
    [RF1,RF2,T1,T2,RD1,RD2]


    % Some BA packages require signed translation. This part solves the
    % sign of the translation. getTsign() computes the sign without the GT
    
    % signM = getTsign(ourR,K_gt,ourT,M);
    scale = sign(scale);
    for i = 1:length(ourT)
        ourT{i} = ourT{i}*scale;
    end

        
    results = struct;
    results.finalTime=finalTime;
    results.optTime=optTime;
    results.RF1 =RF1;
    results.RF2 =RF2;
    results.T1 =T1;
    results.T2 =T2;
    results.RD1 =RD1;
    results.RD2 =RD2;
    dataset



end
