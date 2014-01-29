classdef classDVH2Atlas_MSK_NKI
    properties
        xlsFile_input
        xlsFile_output
        xlsSheet
        xlsRaw % raw data of xls sheet
        flgColInput=true; % patient info. in the input .xls file is listed by columns
        flgColOutput=true; % patient info. in the output .xls file is listed by columns
        
        PatientColName % patient code column name
        PatientCode
        PatientRows % flags of which row has patient code
%         PatientNaN % patients with no DVH data, cannot treat the empty cell as zero in the statistic analysis
        DoseRowName % dose row name
        DoseBins
        DoseCols % flags of columns of dose bins
        ComplicationColName % patient complication code column
        ComplicationGrade % patient complication grade
        ComplicationRows % flags of complication grade
        ComplicationThreshold % threshold of severe complications
        ComplicationFlag % 1 - the complication is severe (grade >= ComplicationThreshold)
        
        lnn = (-1.0:0.1:1)'; % log10(n)
        beta2alpha = 0; % LQ model for dose correction, default is no correction
        FractionNum = 1; % number of fractions, default is single fraction
        GyStep=0.5; % step size for atlas set up
        flgCumulativeDVH = false; % flag for cumulative DVHs, default is differential DVHs
        DVH % DVH data
        EUD % computed EUD
        
        PatientTotal % matrix for total patients
        PatientComp % matrix for patients with complications
        AtlasBins % atlas spacing with step of GyStep
        BetaCumulativeMat % probability that the complication rate is larger than a specific number BetaCumulativeThreshold
        BetaCumulativeThreshold
        BetaInverseMat % The chi-square where the probability is BetaInverseThreshold
        BetaInverseThreshold
        
        LogisticRegressionMatExact % logistic regression model using exact EUD
        LogisticRegressionMatBin % logistic regression model using binned EUD
    end

    methods
        function DAobj = classDVH2Atlas()
        end

        function DAobj=ReadXls(DAobj)
            [~,~,raw] = xlsread(DAobj.xlsFile_input,DAobj.xlsSheet); % read from the xls file
            if ~DAobj.flgColInput % make sure patient info. is listed by columns
                raw=raw';
            end
            DAobj.xlsRaw=raw;
        end
        function DAobj = AssignXlsSheet(DAobj,raw)
            if ~DAobj.flgColInput % make sure patient info. is listed by columns
                raw=raw';
            end
            DAobj.xlsRaw=raw;
        end            
        function DAobj=PatientCodeInXlsraw(DAobj) % extract the column according to column name
            [colcontent,colflg]=ColumnInXlsraw(DAobj.xlsRaw,DAobj.PatientColName); % data in the column
            DAobj.PatientCode=colcontent;
            f=cellfun(@(x) any(isnan(x)), colcontent); % nan cells shall be excluded
            colflg(f)=false; DAobj.PatientRows=colflg; % flags for the data in the column
        end
        function DAobj=ComplicationInXlsraw(DAobj)
            [colcontent,colflg]=ColumnInXlsraw(DAobj.xlsRaw,DAobj.ComplicationColName); % data in the column
            f=cellfun(@(x) any(ischar(x))||any(isnan(x)), colcontent); colflg(f)=false; % nan and character cells shall be excluded
            colflg=colflg&DAobj.PatientRows; % exclude non-patient data
            DAobj.ComplicationRows=colflg; % flags for the data in the column
            DAobj.ComplicationGrade=zeros(size(colcontent)); DAobj.ComplicationGrade(colflg)=cell2mat(colcontent(colflg));
            DAobj.ComplicationFlag=false(size(colcontent)); DAobj.ComplicationFlag(colflg)=DAobj.ComplicationGrade(colflg)>=DAobj.ComplicationThreshold;
        end
        function DAobj=DoseBinInXlsraw(DAobj)
            f=cellfun(@(x) strcmpi(x,DAobj.DoseRowName),DAobj.xlsRaw); % search dose code in the sheet
            [m,~]=find(f); % location of the column
            if length(m)>1
                warning('There are more than one cells containing the same column name, pick the first one.');
            elseif isempty(m)
                error('no column name found, can not continue.');
            end
            
            doses=DAobj.xlsRaw(m(1)+1,:); DAobj.DoseBins=zeros(size(doses));
            f=cellfun(@(x) any(ischar(x))||any(isnan(x)),doses); f=~f; DAobj.DoseBins(f)=cell2mat(doses(f));
            DAobj.DoseCols=f;
        end
        function DAobj=CalculateEUDFromDVH(DAobj)
            % compute EUDs
                % dectect no data patients
                f=cellfun(@(x) any(isnan(x))||any(ischar(x)),DAobj.xlsRaw(:,DAobj.DoseCols));
                [m,~]=find(f); % patients with empty or string DVH data
                DAobj.PatientRows(unique(m))=false; % reomve those patient from EUD computation
%                 DAobj.PatientNaN=false(size(f,1),1); DAobj.PatientNaN(unique(m))=true;

                % prepare differential DVHs
                f=cellfun(@(x) any(isnan(x))||any(ischar(x)),DAobj.xlsRaw);
                f=~f;
                DAobj.DVH=NaN(size(DAobj.xlsRaw)); % the initial of DVH is not zero so that empty cells can be distinguished from cells with zero values
                DAobj.DVH(f)=cell2mat(DAobj.xlsRaw(f));
                if DAobj.flgCumulativeDVH
                    DAobj.DVH(:,1:end-1)=diff(DAobj.DVH,1,2);
                    % check if it is a monotonic decreasing cumulative DVH
                    dvhs=DAobj.DVH(DAobj.PatientRows,DAobj.DoseCols); dvhs(:,end)=0;
                    if any(dvhs(:)>0)
                        error(['The cumulative DVH is not monotonic decreasing. File: ',DAobj.xlsFile_input,' Sheet: ',DAobj.xlsSheet]);
                    end
                    DAobj.DVH=abs(DAobj.DVH); DAobj.DVH(:,~DAobj.DoseCols)=0; % differential volume should be positive, and non-dose columns should be set to zero
                end
                dvhs=DAobj.DVH(:,DAobj.DoseCols);% dvhs(~DAobj.PatientRows,:)=0;

                % transfer doses using LQ model
                dosebins = DAobj.DoseBins(DAobj.DoseCols) .* ( 1 + DAobj.beta2alpha * (DAobj.DoseBins(DAobj.DoseCols) ./ DAobj.FractionNum ) );
                if DAobj.flgCumulativeDVH % dose should be adjusted to be the middle of the bins from the left end
                    dosebins(1:end-1)=dosebins(1:end-1)+diff(dosebins)/2;
                end

                % compute EUD
                nn=10.^DAobj.lnn;
                DAobj.EUD=EUDFromCorrectedDVH( dosebins, dvhs, nn);
        end
        function DAobj=AtlasStatistics(DAobj)
            % complication flags
                if ~isequal(DAobj.PatientRows,DAobj.ComplicationRows)
                    warning('Patients and complication data not consistant');
                end
%                 DAobj.ComplicationFlag=DAobj.ComplicationFlag&DAobj.PatientRows;
%                 DAobj.ComplicationFlag(DAobj.PatientNaN)=false; % confine the flag to available data only
            % for each lnn and every GyStep, compute the total patients and their complications
                DAobj.AtlasBins=(0:DAobj.GyStep:(max(DAobj.EUD(:))+DAobj.GyStep))'; binsnum=length(DAobj.AtlasBins); % DAobj.AtlasBins(1)=DAobj.AtlasBins(1)+eps; % the first element is larger than 0 to exclude no data patients %fix(max(DAobj.EUD(:))/DAobj.GyStep)+1; % total steps(bins) in doses
                DAobj.PatientTotal=zeros(binsnum,length(DAobj.lnn)); DAobj.PatientComp=DAobj.PatientTotal;
                eudall=DAobj.EUD; % eudall(DAobj.PatientNaN,:)=-1; % exclude patients with no data, they definitely can not be treated as with 0 Gy
                eudall(~DAobj.PatientRows,:)=-1; % exclude non-patient rows
                for n=1:length(DAobj.lnn)
                    for m=1:binsnum
                        f=find(eudall(:,n)>=DAobj.AtlasBins(m)); g=find(DAobj.ComplicationFlag(f));
                        DAobj.PatientTotal(m,n)=length(f);
                        DAobj.PatientComp(m,n)=length(g);
                    end
                end
        end
        function DAobj=BetaCumulativeProbability(DAobj)
            DAobj.BetaCumulativeMat=zeros([size(DAobj.PatientTotal) length(DAobj.BetaCumulativeThreshold)]);
            for k=1:length(DAobj.BetaCumulativeThreshold)
                DAobj.BetaCumulativeMat(:,:,k)=betacdf(DAobj.BetaCumulativeThreshold(k),DAobj.PatientComp+1,DAobj.PatientTotal-DAobj.PatientComp+1);
            end
        end
        function DAobj=BetaInverseProbability(DAobj)
            DAobj.BetaInverseMat=zeros([size(DAobj.PatientTotal) length(DAobj.BetaInverseThreshold)]);
            for k=1:length(DAobj.BetaInverseThreshold)
                DAobj.BetaInverseMat(:,:,k)=betainv(DAobj.BetaInverseThreshold(k),DAobj.PatientComp+1,DAobj.PatientTotal-DAobj.PatientComp+1);
            end
        end
        function DAobj=LogisticRegressionAnalysis(DAobj)
            DAobj.LogisticRegressionMatExact = repmat(struct('b',[],'dev',[],'stats',[]),[length(DAobj.lnn),1]);
            DAobj.LogisticRegressionMatBin = repmat(struct('b',[],'dev',[],'stats',[]),[length(DAobj.lnn),1]);
            
            % using exact EUD
            if ~isempty(DAobj.EUD)
                pttotal=ones(size(DAobj.EUD,1),1); pttotal=pttotal(DAobj.PatientRows);
                ptcomp=zeros(size(DAobj.EUD,1),1); ptcomp(DAobj.ComplicationFlag)=1; ptcomp=ptcomp(DAobj.PatientRows);
                for k=1:length(DAobj.lnn)
                    doses=DAobj.EUD(DAobj.PatientRows,k);
                    % regression using exact EUD
                    [b,dev,s]=glmfit(doses,[ptcomp pttotal],'binomial','link','logit');
                    DAobj.LogisticRegressionMatExact(k).b=b; DAobj.LogisticRegressionMatExact(k).dev=dev; DAobj.LogisticRegressionMatExact(k).stats=s;
                end
            end
            % using bins
            dosebins=(DAobj.AtlasBins(1:end-1)+DAobj.AtlasBins(2:end))/2; % dose bins are at the middle of the intervals
            pttotal=ones(DAobj.PatientTotal(1,1),1); % each patient has his own row
            ptcomp=zeros(DAobj.PatientTotal(1,1),1); % allocate space for complication of each patient
            doseval=ptcomp; % allocate space for dose of each patient.
            for k=1:length(DAobj.lnn)
                % the intervals having patients and the num of omplicaitons
                ptintervaltotal=abs(diff(DAobj.PatientTotal(:,k)));
                ptintervalcomp=abs(diff(DAobj.PatientComp(:,k)));
                % generate data set to feed into logistic regression
                n=0; f=find(ptintervaltotal);
                ptcomp(:)=0; doseval(:)=0;
                for m=1:length(f)
                    ptintvnum=ptintervaltotal(f(m)); ptintvcomp=ptintervalcomp(f(m)); % number of patients fall in this interval, and the number of complication within them
                    doseval(n+1:n+ptintvnum)=dosebins(f(m)); % all the patients in the interval share the same dose
                    ptcomp(n+1:n+ptintvcomp)=1;
                    n=n+ptintvnum;
                end
                % logistic regression
                [b,dev,s]=glmfit(doseval,[ptcomp pttotal],'binomial','link','logit');
                DAobj.LogisticRegressionMatBin(k).b=b; DAobj.LogisticRegressionMatBin(k).dev=dev; DAobj.LogisticRegressionMatBin(k).stats=s;
            end
        end
        
        function WriteXls(DAobj)
            warning('off','MATLAB:xlswrite:AddSheet');
            if DAobj.flgColOutput
                % EUD
                xlswrite(DAobj.xlsFile_output,DAobj.lnn',strcat(DAobj.xlsSheet,'_EUD'),'B1');
                if ~isempty(DAobj.PatientCode)
                    xlswrite(DAobj.xlsFile_output,DAobj.PatientCode(DAobj.PatientRows),strcat(DAobj.xlsSheet,'_EUD'),'A2');
                end
                if ~isempty(DAobj.EUD)
                    xlswrite(DAobj.xlsFile_output,DAobj.EUD(DAobj.PatientRows,:),strcat(DAobj.xlsSheet,'_EUD'),'B2');
                end
                
                % total patients at lnn and dose
                xlswrite(DAobj.xlsFile_output,DAobj.lnn', strcat(DAobj.xlsSheet,'_Total'),'B1');
                if ~isempty(DAobj.PatientTotal);
%                     doses = (0:size(DAobj.PatientTotal,1)-1)*DAobj.GyStep;
                    xlswrite(DAobj.xlsFile_output,DAobj.AtlasBins,strcat(DAobj.xlsSheet,'_Total'),'A2');
                    xlswrite(DAobj.xlsFile_output,DAobj.PatientTotal,strcat(DAobj.xlsSheet,'_Total'),'B2');
                end
                
                % complication patients at lnn and dose
                xlswrite(DAobj.xlsFile_output,DAobj.lnn', strcat(DAobj.xlsSheet,'_Comp'),'B1');
                if ~isempty(DAobj.PatientComp)
%                     doses = (0:size(DAobj.PatientComp,1)-1)*DAobj.GyStep;
                    xlswrite(DAobj.xlsFile_output,DAobj.AtlasBins,strcat(DAobj.xlsSheet,'_Comp'),'A2');
                    xlswrite(DAobj.xlsFile_output,DAobj.PatientComp,strcat(DAobj.xlsSheet,'_Comp'),'B2');
                end
                
                % Beta probability
                if ~isempty(DAobj.BetaCumulativeMat)
%                     doses = (0:size(DAobj.PatientComp,1)-1)*DAobj.GyStep;
                    for k=1:length(DAobj.BetaCumulativeThreshold)
                        xlswrite(DAobj.xlsFile_output, DAobj.lnn', strcat(DAobj.xlsSheet,'_prob_',num2str(DAobj.BetaCumulativeThreshold(k))),'B1');
                        xlswrite(DAobj.xlsFile_output,DAobj.AtlasBins,strcat(DAobj.xlsSheet,'_prob_',num2str(DAobj.BetaCumulativeThreshold(k))),'A2');
                        xlswrite(DAobj.xlsFile_output,DAobj.BetaCumulativeMat(:,:,k),strcat(DAobj.xlsSheet,'_prob_',num2str(DAobj.BetaCumulativeThreshold(k))),'B2');
                    end
                end
                if ~isempty(DAobj.BetaInverseMat)
%                     doses = (0:size(DAobj.PatientComp,1)-1)*DAobj.GyStep;
                    for k=1:length(DAobj.BetaInverseThreshold)
                        xlswrite(DAobj.xlsFile_output,DAobj.lnn', strcat(DAobj.xlsSheet,'_Low_',num2str(DAobj.BetaInverseThreshold(k))),'B1');
                        xlswrite(DAobj.xlsFile_output,DAobj.AtlasBins,strcat(DAobj.xlsSheet,'_Low_',num2str(DAobj.BetaInverseThreshold(k))),'A2');
                        xlswrite(DAobj.xlsFile_output,DAobj.BetaCumulativeMat(:,:,k),strcat(DAobj.xlsSheet,'_Low_',num2str(DAobj.BetaInverseThreshold(k))),'B2');
                    end
                end
            else
                % EUD
                xlswrite(DAobj.xlsFile_output, DAobj.lnn, strcat(DAobj.xlsSheet,'_EUD'),'A2');
                if ~isempty(DAobj.PatientCode)
                    xlswrite(DAobj.xlsFile_output,DAobj.PatientCode(DAobj.PatientRows)',strcat(DAobj.xlsSheet,'_EUD'),'B1');
                end
                if ~isempty(DAobj.EUD)
                    xlswrite(DAobj.xlsFile_output,DAobj.EUD(DAobj.PatientRows,:)',strcat(DAobj.xlsSheet,'_EUD'),'B2');
                end
                
                % total patients at lnn and dose
                xlswrite(DAobj.xlsFile_output, DAobj.lnn, strcat(DAobj.xlsSheet,'_Total'),'A2');
                if ~isempty(DAobj.PatientTotal);
%                     doses = (0:size(DAobj.PatientTotal,1)-1)*DAobj.GyStep;
                    xlswrite(DAobj.xlsFile_output,DAobj.AtlasBins',strcat(DAobj.xlsSheet,'_Total'),'B1');
                    xlswrite(DAobj.xlsFile_output,DAobj.PatientTotal',strcat(DAobj.xlsSheet,'_Total'),'B2');
                end
                
                % complication patients at lnn and dose
                xlswrite(DAobj.xlsFile_output, DAobj.lnn, strcat(DAobj.xlsSheet,'_Comp'),'A2');
                if ~isempty(DAobj.PatientComp)
%                     doses = (0:size(DAobj.PatientComp,1)-1)*DAobj.GyStep;
                    xlswrite(DAobj.xlsFile_output,DAobj.AtlasBins',strcat(DAobj.xlsSheet,'_Comp'),'B1');
                    xlswrite(DAobj.xlsFile_output,DAobj.PatientComp',strcat(DAobj.xlsSheet,'_Comp'),'B2');
                end
                
                % Beta probability
                if ~isempty(DAobj.BetaCumulativeMat)
%                     doses = (0:size(DAobj.PatientComp,1)-1)*DAobj.GyStep;
                    for k=1:length(DAobj.BetaCumulativeThreshold)
                        xlswrite(DAobj.xlsFile_output, DAobj.lnn, strcat(DAobj.xlsSheet,'_prob_',num2str(DAobj.BetaCumulativeThreshold(k))),'A2');
                        xlswrite(DAobj.xlsFile_output,DAobj.AtlasBins',strcat(DAobj.xlsSheet,'_prob_',num2str(DAobj.BetaCumulativeThreshold(k))),'B1');
                        xlswrite(DAobj.xlsFile_output,DAobj.BetaCumulativeMat(:,:,k)',strcat(DAobj.xlsSheet,'_prob_',num2str(DAobj.BetaCumulativeThreshold(k))),'B2');
                    end
                end
                if ~isempty(DAobj.BetaInverseMat)
%                     doses = (0:size(DAobj.PatientComp,1)-1)*DAobj.GyStep;
                    for k=1:length(DAobj.BetaInverseThreshold)
                        xlswrite(DAobj.xlsFile_output, DAobj.lnn, strcat(DAobj.xlsSheet,'_Low_',num2str(DAobj.BetaInverseThreshold(k))),'A2');
                        xlswrite(DAobj.xlsFile_output,DAobj.AtlasBins',strcat(DAobj.xlsSheet,'_Low_',num2str(DAobj.BetaInverseThreshold(k))),'B1');
                        xlswrite(DAobj.xlsFile_output,DAobj.BetaInverseMat(:,:,k)',strcat(DAobj.xlsSheet,'_Low_',num2str(DAobj.BetaInverseThreshold(k))),'B2');
                    end
                end
            end
            warning('on','MATLAB:xlswrite:AddSheet');
        end
        function AtlasFig(DAobj)
            if isempty(DAobj.PatientTotal)
                disp('Empty member "PatientTotal", cannot display its figure.'); return;
            end
            doses=(0:size(DAobj.PatientTotal,1)-1)*DAobj.GyStep; x=mod(doses,5)==0;
            [xx,yy]=ndgrid(doses(x),1:length(DAobj.lnn)); xx=num2cell(xx); yy=num2cell(yy);
            str1=(DAobj.PatientComp(x,:)); str2=(DAobj.PatientTotal(x,:)); str=arrayfun(@(a,b) strcat(num2str(a),'/',num2str(b)),str1,str2,'UniformOutpu',false);
            cellfun(@(a,b,c) text(a,b,c,'fontsize',16),xx,yy,str);
            set(gca,'XLim',[xx{1,1}-1,xx{end,1}+2]); set(gca,'YLim',[1-1,length(DAobj.lnn)+1]);
            set(gca,'YTick',1:length(DAobj.lnn)); set(gca,'YTickLabel',DAobj.lnn);
            pos=get(gcf,'Position'); set(gcf,'Position',[pos(1:2), 900/3*4, 900]);
            title([DAobj.xlsSheet,', number of patients and numbers with severe pneumonitis treated to at least a given EUD'],'fontsize',16);
            xlabel('EUD doses (Gy)','fontsize',16); ylabel('log n','fontsize',16); set(gca,'fontsize',16);
        end
        function ProbabilityFig(DAobj)
            if isempty(DAobj.BetaCumulativeMat)
                disp('Empty member (BetaCumulativeMat), can not display its figure.'); return;
            end
            img=1-DAobj.BetaCumulativeMat'; img(1,end)=1; img(2,end)=0;
            contourf(img); %contourf(flipud(rot90(DAobj.BetaCumulativeMat,1)));
            xsteps=0:DAobj.GyStep:max(DAobj.EUD(:)); xtickvec=10:10:xsteps(end); xtickstep=(size(DAobj.BetaCumulativeMat,1)/xsteps(end)*10);
            set(gca,'XTick',xtickstep:xtickstep:size(DAobj.BetaCumulativeMat,1)); set(gca,'XTickLabel',xtickvec);
            ytickvec=DAobj.lnn; ytickstep=size(DAobj.BetaCumulativeMat,2)/length(ytickvec);
            set(gca,'YTick',ytickstep:ytickstep:size(DAobj.BetaCumulativeMat,2)); set(gca,'YTickLabel',ytickvec);
            %             pmin=min(DAobj.BetaCumulativeMat(:)); pmax=max(DAobj.BetaCumulativeMat(:)); %bstep=(pmax-pmin)/10*256;
            colorbar; %colorbar('YTickLabel',round([pmin:(pmax-pmin)/9:pmax]*10)/10);
            title([DAobj.xlsSheet,', the probability that observed complications arise from true rate > 20%']);
            xlabel('EUD doses (Gy)'); ylabel('log n');
        end
        function LowProbabilityFig(DAobj)
            if isempty(DAobj.BetaInverseMat)
                disp('Empty member (BetaInverseMat), can not display its figure.'); return;
            end
            img=DAobj.BetaInverseMat'; img(1,end)=1; img(2,end)=0;
            contourf(img); %contourf(flipud(rot90(DAobj.BetaInverseMat,1)));
            xsteps=0:DAobj.GyStep:max(DAobj.EUD(:)); xtickvec=10:10:xsteps(end); xtickstep=(size(DAobj.BetaInverseMat,1)/xsteps(end)*10);
            set(gca,'XTick',xtickstep:xtickstep:size(DAobj.BetaInverseMat,1)); set(gca,'XTickLabel',xtickvec);
            ytickvec=DAobj.lnn; ytickstep=size(DAobj.BetaInverseMat,2)/length(ytickvec);
            set(gca,'YTick',ytickstep:ytickstep:size(DAobj.BetaInverseMat,2)); set(gca,'YTickLabel',ytickvec);
            colorbar;
            title([DAobj.xlsSheet,', lower 68% confidence limit on complication probability']);
            xlabel('EUD doses (Gy)'); ylabel('log n');
        end
        
        function DAobj=CombineAtlas(DAobj,DAobj1,DAobj2)
            if ~isequal(DAobj1.GyStep,DAobj2.GyStep)
                error('Bin steps in the Atlas in the two DVH2Atlas objects are different, can not merge the atlas');
            end
            DAobj=DAobj1;
            DAobj.PatientCode=[DAobj.PatientCode; DAobj2.PatientCode];
            DAobj.PatientRows=[DAobj.PatientRows; DAobj2.PatientRows];
%             DAobj.PatientNaN=[DAobj.PatientNaN; DAobj2.PatientNaN];
            DAobj.ComplicationGrade=[DAobj.ComplicationGrade;DAobj2.ComplicationGrade];
            DAobj.ComplicationRows=[DAobj.ComplicationRows;DAobj2.ComplicationRows];
            DAobj.ComplicationFlag=[DAobj.ComplicationFlag; DAobj2.ComplicationFlag];
            
            if isempty(DAobj.EUD) || isempty(DAobj2.EUD) || sum(abs(DAobj.lnn-DAobj2.lnn))>1.0e-9
                DAobj.EUD=[];
            else
                DAobj.EUD=[DAobj.EUD;DAobj2.EUD];
            end
            
            DAobj.AtlasBins=unique([DAobj.AtlasBins;DAobj2.AtlasBins]);
            sizemsk=size(DAobj.PatientTotal); sizenki=size(DAobj2.PatientTotal);
            DAobj.PatientTotal=zeros(max(sizemsk,sizenki));
            DAobj.PatientTotal(1:sizemsk(1),1:sizemsk(2))=DAobj1.PatientTotal;
            DAobj.PatientTotal(1:sizenki(1),1:sizenki(2))=DAobj.PatientTotal(1:sizenki(1),1:sizenki(2))+DAobj2.PatientTotal;
            DAobj.PatientComp=zeros(max(sizemsk,sizenki));
            DAobj.PatientComp(1:sizemsk(1),1:sizemsk(2))=DAobj1.PatientComp;
            DAobj.PatientComp(1:sizenki(1),1:sizenki(2))=DAobj.PatientComp(1:sizenki(1),1:sizenki(2))+DAobj2.PatientComp;
            
            % compute probabilities
            DAobj=DAobj.BetaCumulativeProbability();
            DAobj=DAobj.BetaInverseProbability();
            DAobj=DAobj.LogisticRegressionAnalysis();
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
