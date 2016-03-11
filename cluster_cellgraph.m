% folder names
foldStr = {'Resistant','Susceptible','Supersusceptible'};

% assuming value of number of regions within an image
regions = 4;
medEcc = 0.5;

% create bash file to run the cell graph code
fp = fopen('bash_3classes.sh','w');

% loop through MAT files in each folder to perform analysis
for ind = 1 : size(foldStr,2)
    
    %  find all .MAT files within this folder
    mats = dir([foldStr{ind}, '/*.mat']);
    
    for iter = 1 : size(mats,1)
        % load .MAT file
        feats = load([foldStr{ind}, '/' mats(iter).name]);
        temp = fieldnames(feats);
        feats = feats.(temp{1});
        feats = feats(:,[1:8, 10, 11]);
        
        % histogram the red color channel into specified number of regions
        [~,~,binRed] = histcounts(feats(:,4),regions);
        % hist(feats(:,6),regions);
        
        % split each region based on size and elongation into 2 different
        % sub populations
        for ind2 = 1 : regions
            cellsRegion = feats(binRed == ind2,:);
            
            % split based on size
            medSize = median(cellsRegion(:,1));
            smallSz = cellsRegion(cellsRegion(:,1) < medSize,:);
            meanDistSm = mean(smallSz(:,10));
            largeSz = cellsRegion(cellsRegion(:,1) >= medSize,:);
            meanDistLg = mean(largeSz(:,10));
            
            % save txt files and add entries to bash file
            fileName = [foldStr{ind} '/cells_small_' mats(iter).name '_' ...
                num2str(ind2) '_region.txt'];
            fileid = fopen(fileName,'w');           
            fprintf(fileid,'%f %f\n',smallSz(:,[2 3]));            
            fclose(fileid);
            
            fprintf(fp,'./GraphMetrics %s graphs_small_%s_%d.txt metrics_small_%s_%d.txt %d\n', ...
                fileName, mats(iter).name, ind2, mats(iter).name, ind2, meanDistSm);
            
            fileName = [foldStr{ind} '/cells_large_' mats(iter).name '_' ...
                num2str(ind2) '_region.txt'];
            fileid = fopen(fileName,'w');           
            fprintf(fileid,'%f %f\n',largeSz(:,[2 3]));            
            fclose(fileid);
            
            fprintf(fp,'./GraphMetrics %s graphs_large_%s_%d.txt metrics_large_%s_%d.txt %d\n', ...
                fileName, mats(iter).name, ind2, mats(iter).name, ind2, meanDistLg);
            
            % split based on eccentricity
            smallEcc = cellsRegion(cellsRegion(:,9) < medEcc,:);
            meanEccSm = mean(smallEcc(:,10));
            largeEcc = cellsRegion(cellsRegion(:,9) >= medEcc,:);
            meanEccLg = mean(largeEcc(:,10));
            
            % save txt files and add entries to bash file
            fileName = [foldStr{ind} '/cells_lessecc_' mats(iter).name '_' ...
                num2str(ind2) '_region.txt'];
            fileid = fopen(fileName,'w');           
            fprintf(fileid,'%f %f\n',smallEcc(:,[2 3]));            
            fclose(fileid);
            
            fprintf(fp,'./GraphMetrics %s graphs_lowecc_%s_%d.txt metrics_lowecc_%s_%d.txt %d\n', ...
                fileName, mats(iter).name, ind2, mats(iter).name, ind2, meanEccSm);
            
            fileName = [foldStr{ind} '/cells_moreecc_' mats(iter).name '_' ...
                num2str(ind2) '_region.txt'];
            fileid = fopen(fileName,'w');           
            fprintf(fileid,'%f %f\n',largeEcc(:,[2 3]));            
            fclose(fileid);
            
            fprintf(fp,'./GraphMetrics %s graphs_hiecc_%s_%d.txt metrics_hiecc_%s_%d.txt %d\n', ...
                fileName, mats(iter).name, ind2, mats(iter).name, ind2, meanEccLg);
        end
        
        % nothing to be done for green region, maybe reconsider later (?)
        % [NG,~,binGreen] = histcounts(feats(:,5),regions);
    end    
end
fclose(fp); 

% for i = 1 : K
%     fileName = ['cells_' num2str(i) '_region.txt'];
%     fileid = fopen(fileName,'w');
%     fprintf(fileid,'%f %f %f\n',cellsRegion);
%     fclose(fileid);
%     str = ['../../GraphMetrics ' fileName ' graphs_' num2str(i) '.txt' ...
% 
%         ' metrics_' num2str(i) '.txt ' num2str(meanRegion)];
%     unix(str);
% end