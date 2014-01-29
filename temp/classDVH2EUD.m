classdef DVH2EUD
    properties
        xlsFile
        xlsSheet
        xlsDESC
        xlsRaw % raw data of xls sheet
        flgColInput=true; % patient info. in the input .xls file is listed by columns
        flgColOutput=true; % patient info. in the output .xls file is listed by columns
        
        rows % vector showing the rows of DVH info in the sheet. Fist row is dose info, each other row is DVH info of a patient
        cols % vector showing the columns (positions) of DVH info
        
        PatientColName % patient code column name
        PatientCode
        PatientRows
        PatientNaN % patients with no DVH data, cannot treat the empty cell as zero in the statistic analysis
        DoseRowName % dose row name
        DoseBins
        DoseCols
        ComplicationColName % patient complication code column
        ComplicationGrade % patient complication grade
        ComplicationRows
        ComplicationThreshold % threshold of severe complications
        
        lnn = (-1.0:0.1:1)';
        nn = 10.^(-1.0:0.1:1)';
        beta2alpha = 0;
        FractionNum = 1;
        GyStep=0.5;
        flgCumulativeDVH = false;
        DVH
        EUD
        
        PatientTotal % matrix for total patients
        PatientComp % matrix for patients with complications
        BetaCumulativeMat
        BetaCumulativeThreshold
        BetaInverseMat
        BetaInverseThreshold
    end

    methods
        function DEobj = DVH2EUD()
        end
        function DEobj=set.lnn(DEobj,lnn)
            DEobj.lnn=lnn(:);
            DEobj.nn=10.^DEobj.lnn;
        end
        function DEobj=set.xlsFile(DEobj,xlsFile)
            [xlsinfo1,xlsinfo2]=xlsfinfo(xlsFile);
            if ~strcmp('Microsoft Excel Spreadsheet',xlsinfo1)
                error('not a .xls file');
            end
            DEobj.xlsFile=xlsFile;
            DEobj.xlsDESC=xlsinfo2;
        end
        function DEobj=set.xlsSheet(DEobj,xlsSheet)
            if ~any(strcmp(DEobj.xlsDESC,xlsSheet))
                error(['not a correct sheet for the file ',DEobj.xlsFile]);
            end
            DEobj.xlsSheet=xlsSheet;
            % read the sheet
            [~,~,raw] = xlsread(DEobj.xlsFile,xlsSheet); % read from the xls file
            if ~DEobj.flgColInput % make sure patient info. is listed by columns
                raw=raw';
            end
            DEobj.xlsRaw=raw;
        end
        function DEobj=set.flgColInput(DEobj,flg)
            if DEobj.flgColInput==flg
                return;
            end
            DEobj.flgColInput=flg;
            DEobj.xlsRaw=DEobj.xlsRaw';
        end

        function DEobj=PatientCodeInXlsraw(DEobj) % extract the column according to column name
            [colcontent,colflg]=ColumnInXlsraw(DEobj.xlsRaw,DEobj.PatientColName); % data in the column
            DEobj.PatientCode=colcontent;
            f=cellfun(@(x) any(isnan(x)), colcontent); % nan cells shall be excluded
            colflg(f)=false; DEobj.PatientRows=colflg; % flags for the data in the column
        end
        function DEobj=ComplicationInXlsraw(DEobj)
            [colcontent,colflg]=ColumnInXlsraw(DEobj.xlsRaw,DEobj.ComplicationColName); % data in the column
            f=cellfun(@(x) ischar(x)||any(isnan(x)), colcontent); % nan and character cells shall be excluded
            colflg(f)=false; DEobj.ComplicationRows=colflg; % flags for the data in the column
            DEobj.ComplicationGrade=zeros(size(colcontent)); DEobj.ComplicationGrade(colflg)=cell2mat(colcontent(colflg));
        end
        function DEobj=DoseBinInXlsraw(DEobj)
            f=cellfun(@(x) strcmpi(x,DEobj.DoseRowName),DEobj.xlsRaw); % search dose code in the sheet
            [m,~]=find(f); % location of the column
            if length(m)>1
                warning('There are more than one cells containing the same column name, pick the first one.');
            elseif isempty(m)
                error('no column name found, can not continue.');
            end
            
            doses=DEobj.xlsRaw(m(1)+1,:); DEobj.DoseBins=zeros(size(doses));
            f=cellfun(@(x) any(ischar(x))|any(isnan(x)),doses); f=~f; DEobj.DoseBins(f)=cell2mat(doses(f));
            DEobj.DoseCols=f;
        end
        function DEobj=CalculateEUDFromDVH(DEobj)
%             % adjust data according to the flgColInput
%                 if DEobj.flgColInput
%                     dataRows=DEobj.rows; dataCols=DEobj.cols; xlsdata=DEobj.xlsRaw;
%                 else
%                     dataRows=DEobj.cols; dataCols=DEobj.rows; xlsdata=DEobj.xlsRaw';
%                 end
%             % retrieve patient code
%                 flg=cellfun(@(x) strcmpi(x,DEobj.PatienttCol),DEobj.xlsRaw);
%                 if sum(flg)>1
%                     disp('There are more than one cells containing patient code column name, pick the first one to retrieve patient code.');
%                 end
%                 if sum(flg)==0
%                     warning('no patient code found, can not continue.');
%                     return;
%                 end
%                 [m,n]=find(flg); % the location of the cell of patient code
%                 DEobj.PatientCode=DEobj.xlsRaw(m(1),dataCols(2:end));

            % compute EUDs
                % dectect no data patients
                f=cellfun(@(x) any(isnan(x))||any(ischar(x)),DEobj.xlsRaw(:,DEobj.DoseCols));
                [m,~]=find(f); % patients with no DVH data
                DEobj.PatientNaN=f(:,1); DEobj.PatientNaN(:)=false; DEobj.PatientNaN(unique(m))=true;
                % prepare differential DVHs
                f=cellfun(@(x) any(isnan(x))||any(ischar(x)),DEobj.xlsRaw);
                f=~f;
                DEobj.DVH=zeros(size(DEobj.xlsRaw));
                DEobj.DVH(f)=cell2mat(DEobj.xlsRaw(f));
%                 [m,~]=find(f); m=unique(m);
%                 f=find(DEobj.PatientRows); DEobj.PatientRows(f(m))=false;
                % compute EUDs
%                 DEobj.DVH = cell2mat( DEobj.xlsRaw(DEobj.PatientRows,DEobj.DoseCols) );
                if DEobj.flgCumulativeDVH
                    DEobj.DVH(:,1:end-1)=diff(DEobj.DVH,1,2);
                    % check if it is cumulative DVH
                    dvhs=DEobj.DVH(DEobj.PatientRows,DEobj.DoseCols); dvhs(:,end)=0;
                    if any(dvhs(:)>0)
                        error(['The DVH is not cumulative in file ',DEobj.xlsFile,' sheet: ',DEobj.xlsSheet,'\n but was set as']);
                    end
                    DEobj.DVH=abs(DEobj.DVH); DEobj.DVH(:,~DEobj.DoseCols)=0; % differential volume should be positive, and non-dose columns should be set to zero
                end
                dvhs=DEobj.DVH(:,DEobj.DoseCols); dvhs(~DEobj.PatientRows,:)=0;
                % transfer doses using LQ model
                dosebins = DEobj.DoseBins(DEobj.DoseCols) .* ( 1 + DEobj.beta2alpha * (DEobj.DoseBins(DEobj.DoseCols) ./ DEobj.FractionNum ) );
                % compute EUD
                DEobj.EUD=EUDFromCorrectedDVH( dosebins, dvhs, DEobj.nn);
%                 celldose = DEobj.xlsRaw( dataRows, dataCols(1) ); 
%                 f=cellfun(@(x) ischar(x),celldose); DEobj.Doses=zeros(size(celldose)); DEobj.Doses(~f)=cell2mat(celldose(~f));
%                 celldose = DEobj.xlsRaw( dataRows, dataCols(2:end) );
%                 f=cellfun(@(x) ischar(x), celldose); DEobj.DVH=zeros(size(celldose)); DEobj.DVH(~f)=cell2mat(celldose(~f));
%                 DEobj.EUD = EUDFromCorrectedDVH( DEobj.Doses, DEobj.DVH, DEobj.nn);
        end
        function DEobj=AtlasStatistics(DEobj)
            % prepare
                DEobj.PatientTotal=[]; DEobj.PatientComp=[];
            % complication flags
                if ~isequal(DEobj.PatientRows,DEobj.ComplicationRows)
                    warning('Patients and complication data not consistant');
                end
                flgcom=DEobj.ComplicationGrade>=DEobj.ComplicationThreshold; flgcom=flgcom&DEobj.ComplicationRows&~DEobj.PatientNaN&DEobj.PatientRows; % confine the flag to available data only
            % for each lnn and every GyStep, compute the total patients and their complications
                stps=0:DEobj.GyStep:max(DEobj.EUD(:)); stpstotal=length(stps); % stps(1)=stps(1)+eps; % the first element is larger than 0 to exclude no data patients %fix(max(DEobj.EUD(:))/DEobj.GyStep)+1; % total steps(bins) in doses
                DEobj.PatientTotal=zeros(stpstotal,length(DEobj.lnn)); DEobj.PatientComp=DEobj.PatientTotal;
                eudall=DEobj.EUD; % eudall(DEobj.PatientNaN,:)=-1; % exclude patients with no data, they definitely can not be treated as with 0 Gy
                % eudall(~DEobj.PatientRows,:)=-1; % exclude non-patient rows
                for n=1:length(DEobj.lnn)
                    for m=1:stpstotal
                        f=find(eudall(:,n)>=stps(m)); g=find(flgcom(f));
                        DEobj.PatientTotal(m,n)=length(f);
                        DEobj.PatientComp(m,n)=length(g);
                    end
                end
        end
        function DEobj=BetaCumulativeProbability(DEobj)
            DEobj.BetaCumulativeMat=zeros([size(DEobj.PatientTotal) length(DEobj.BetaCumulativeThreshold)]);
            for k=1:length(DEobj.BetaCumulativeThreshold)
                DEobj.BetaCumulativeMat(:,:,k)=betacdf(DEobj.BetaCumulativeThreshold(k),DEobj.PatientComp+1,DEobj.PatientTotal-DEobj.PatientComp+1);
            end
        end
        function DEobj=BetaInverseProbability(DEobj)
            DEobj.BetaInverseMat=zeros([size(DEobj.PatientTotal) length(DEobj.BetaInverseThreshold)]);
            for k=1:length(DEobj.BetaInverseThreshold)
                DEobj.BetaInverseMat(:,:,k)=betainv(DEobj.BetaInverseThreshold(k),DEobj.PatientComp+1,DEobj.PatientTotal-DEobj.PatientComp+1);
            end
        end
        function xlsWrite(DEobj)
            warning('off','MATLAB:xlswrite:AddSheet');
            if DEobj.flgColOutput
                % EUD
                xlswrite(strcat(DEobj.xlsFile,'_tom'),[DEobj.lnn'; DEobj.nn'],strcat(DEobj.xlsSheet,'_EUD'),'B1');
                if ~isempty(DEobj.PatientCode)
                    xlswrite(strcat(DEobj.xlsFile,'_tom'),DEobj.PatientCode(DEobj.PatientRows),strcat(DEobj.xlsSheet,'_EUD'),'A3');
                end
                if ~isempty(DEobj.EUD)
                    xlswrite(strcat(DEobj.xlsFile,'_tom'),DEobj.EUD(DEobj.PatientRows,:),strcat(DEobj.xlsSheet,'_EUD'),'B3');
                end
                
                % total patients at lnn and dose
                xlswrite(strcat(DEobj.xlsFile,'_tom'),[DEobj.lnn'; DEobj.nn'],strcat(DEobj.xlsSheet,'_Total'),'B1');
                if ~isempty(DEobj.PatientTotal);
                    doses = (0:size(DEobj.PatientTotal,1)-1)*DEobj.GyStep;
                    xlswrite(strcat(DEobj.xlsFile,'_tom'),doses',strcat(DEobj.xlsSheet,'_Total'),'A3');
                    xlswrite(strcat(DEobj.xlsFile,'_tom'),DEobj.PatientTotal,strcat(DEobj.xlsSheet,'_Total'),'B3');
                end
                
                % complication patients at lnn and dose
                xlswrite(strcat(DEobj.xlsFile,'_tom'),[DEobj.lnn'; DEobj.nn'],strcat(DEobj.xlsSheet,'_Comp'),'B1');
                if ~isempty(DEobj.PatientComp)
                    doses = (0:size(DEobj.PatientComp,1)-1)*DEobj.GyStep;
                    xlswrite(strcat(DEobj.xlsFile,'_tom'),doses',strcat(DEobj.xlsSheet,'_Comp'),'A3');
                    xlswrite(strcat(DEobj.xlsFile,'_tom'),DEobj.PatientComp,strcat(DEobj.xlsSheet,'_Comp'),'B3');
                end
                
                % Beta probability
                if ~isempty(DEobj.BetaCumulativeMat)
                    doses = (0:size(DEobj.PatientComp,1)-1)*DEobj.GyStep;
                    for k=1:length(DEobj.BetaCumulativeThreshold)
                        xlswrite(strcat(DEobj.xlsFile,'_tom'),[DEobj.lnn, DEobj.nn],strcat(DEobj.xlsSheet,'_prob_',num2str(DEobj.BetaCumulativeThreshold(k))),'B1');
                        xlswrite(strcat(DEobj.xlsFile,'_tom'),doses',strcat(DEobj.xlsSheet,'_prob_',num2str(DEobj.BetaCumulativeThreshold(k))),'A3');
                        xlswrite(strcat(DEobj.xlsFile,'_tom'),DEobj.BetaCumulativeMat(:,:,k),strcat(DEobj.xlsSheet,'_prob_',num2str(DEobj.BetaCumulativeThreshold(k))),'B3');
                    end
                end
                if ~isempty(DEobj.BetaInverseMat)
                    doses = (0:size(DEobj.PatientComp,1)-1)*DEobj.GyStep;
                    for k=1:length(DEobj.BetaInverseThreshold)
                        xlswrite(strcat(DEobj.xlsFile,'_tom'),[DEobj.lnn, DEobj.nn],strcat(DEobj.xlsSheet,'_Low_',num2str(DEobj.BetaInverseThreshold(k))),'B1');
                        xlswrite(strcat(DEobj.xlsFile,'_tom'),doses',strcat(DEobj.xlsSheet,'_Low_',num2str(DEobj.BetaInverseThreshold(k))),'A3');
                        xlswrite(strcat(DEobj.xlsFile,'_tom'),DEobj.BetaCumulativeMat(:,:,k),strcat(DEobj.xlsSheet,'_Low_',num2str(DEobj.BetaInverseThreshold(k))),'B3');
                    end
                end
            else
                % EUD
                xlswrite(strcat(DEobj.xlsFile,'_tom'),[DEobj.lnn, DEobj.nn],strcat(DEobj.xlsSheet,'_EUD'),'A2');
                if ~isempty(DEobj.PatientCode)
                    xlswrite(strcat(DEobj.xlsFile,'_tom'),DEobj.PatientCode(DEobj.PatientRows)',strcat(DEobj.xlsSheet,'_EUD'),'C1');
                end
                if ~isempty(DEobj.EUD)
                    xlswrite(strcat(DEobj.xlsFile,'_tom'),DEobj.EUD(DEobj.PatientRows,:)',strcat(DEobj.xlsSheet,'_EUD'),'C2');
                end
                
                % total patients at lnn and dose
                xlswrite(strcat(DEobj.xlsFile,'_tom'),[DEobj.lnn, DEobj.nn],strcat(DEobj.xlsSheet,'_Total'),'A2');
                if ~isempty(DEobj.PatientTotal);
                    doses = (0:size(DEobj.PatientTotal,1)-1)*DEobj.GyStep;
                    xlswrite(strcat(DEobj.xlsFile,'_tom'),doses,strcat(DEobj.xlsSheet,'_Total'),'C1');
                    xlswrite(strcat(DEobj.xlsFile,'_tom'),DEobj.PatientTotal',strcat(DEobj.xlsSheet,'_Total'),'C2');
                end
                
                % complication patients at lnn and dose
                xlswrite(strcat(DEobj.xlsFile,'_tom'),[DEobj.lnn, DEobj.nn],strcat(DEobj.xlsSheet,'_Comp'),'A2');
                if ~isempty(DEobj.PatientComp)
                    doses = (0:size(DEobj.PatientComp,1)-1)*DEobj.GyStep;
                    xlswrite(strcat(DEobj.xlsFile,'_tom'),doses,strcat(DEobj.xlsSheet,'_Comp'),'C1');
                    xlswrite(strcat(DEobj.xlsFile,'_tom'),DEobj.PatientComp',strcat(DEobj.xlsSheet,'_Comp'),'C2');
                end
                
                % Beta probability
                if ~isempty(DEobj.BetaCumulativeMat)
                    doses = (0:size(DEobj.PatientComp,1)-1)*DEobj.GyStep;
                    for k=1:length(DEobj.BetaCumulativeThreshold)
                        xlswrite(strcat(DEobj.xlsFile,'_tom'),[DEobj.lnn, DEobj.nn],strcat(DEobj.xlsSheet,'_prob_',num2str(DEobj.BetaCumulativeThreshold(k))),'A2');
                        xlswrite(strcat(DEobj.xlsFile,'_tom'),doses,strcat(DEobj.xlsSheet,'_prob_',num2str(DEobj.BetaCumulativeThreshold(k))),'C1');
                        xlswrite(strcat(DEobj.xlsFile,'_tom'),DEobj.BetaCumulativeMat(:,:,k)',strcat(DEobj.xlsSheet,'_prob_',num2str(DEobj.BetaCumulativeThreshold(k))),'C2');
                    end
                end
                if ~isempty(DEobj.BetaInverseMat)
                    doses = (0:size(DEobj.PatientComp,1)-1)*DEobj.GyStep;
                    for k=1:length(DEobj.BetaInverseThreshold)
                        xlswrite(strcat(DEobj.xlsFile,'_tom'),[DEobj.lnn, DEobj.nn],strcat(DEobj.xlsSheet,'_Low_',num2str(DEobj.BetaInverseThreshold(k))),'A2');
                        xlswrite(strcat(DEobj.xlsFile,'_tom'),doses,strcat(DEobj.xlsSheet,'_Low_',num2str(DEobj.BetaInverseThreshold(k))),'C1');
                        xlswrite(strcat(DEobj.xlsFile,'_tom'),DEobj.BetaInverseMat(:,:,k)',strcat(DEobj.xlsSheet,'_Low_',num2str(DEobj.BetaInverseThreshold(k))),'C2');
                    end
                end
            end
            warning('on','MATLAB:xlswrite:AddSheet');
        end
        
    end
end


function [colcontent,colflg]=ColumnInXlsraw(xlsraw,colname)
% extract the column starting from the row under column name
    f=cellfun(@(x) strcmpi(x,colname), xlsraw);
    [m,n]=find(f); % location of the column
    if length(m)>1
        warning('There are more than one cells containing the same column name, pick the first one.');
    elseif isempty(m)
        error('no column name found, can not continue.');
    end
    colcontent = xlsraw(:,n); colflg=true(size(xlsraw,1),1); colflg(1:m(1))=false;
end


function eud=EUDFromCorrectedDVH(doserow,volumemat,n)
% doserow -- one row of dose info
% volumemat -- rows coresponds to patients' volume info
% EUD computation according to (log10 n), which changes from -1 to 1 with step 0.1

% prepare
    [dimx,dimy]=size(volumemat);

% partial volume computation
    v=sum(volumemat,2);
    if any(v~=1) % not partial volume, transfer it
        volumemat=volumemat./repmat(v,[1,dimy]);
    end
    
% Compute EUD for different ns
    % sum dose under each n to form EUD
    eud=zeros(dimx,length(n));
    for k=1:length(n)
        d1=doserow.^(1/n(k)); d2=repmat(d1,[dimx,1]); % (di)^(1/n)
        d2=d2.*volumemat; % (di)^(1/n) * vi
        eud(:,k)=(sum(d2,2)).^n(k);
    end
end
