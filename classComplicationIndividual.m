classdef classComplicationIndividual
    properties
        PatientID
        Gender
        BirthDate % used to compute the patient age at some time point
        DeathDate % used to compute the survival time
        
        DosePrescription
        DoseBins_org
        VolDiff
        VolCum
        Beta2Alpha
        DoseBins_LQ
        FxNum
        
%         lnn = (-1:0.1:1)'; % log10(n) for EUD computation
        lgn % log10(n)
        nLyman % parameter n in Lyman model
        EUD
        
%         DoseStep
%         VolStep
        
        Stage % stage before the treatment
        BaselineDate % the date of the start point (in number)
        CompOccurDate % the date of complication occurance (in number)
        CompGrade % patient complication grade
        CompThreshold % threshold of severe complications
        flgCensor % flag for censor, 1 - censored
        LastFollowupDate % the date of the last followup
        RelapseDate % the date of relapse (in number)
        
    end
    
    methods % operations
        function CIobj = classComplicationIndividual()
        end
        function CIobj = set.PatientID(CIobj,ptid)
            CIobj.PatientID = num2str(ptid);
        end
        function CIobj = set.lgn(CIobj,lgn)
            CIobj.lgn = lgn(:);
        end
        function CIobj = set.nLyman(CIobj,n)
            CIobj.nLyman = n(:);
        end
        function CIobj = set.DoseBins_org(CIobj,dosebins)
            if any(isnan(dosebins))
                error('there is nan element when atempting to assigning member DoseBins_org, in classComplicationIndividual');
            end
            CIobj.DoseBins_org = dosebins(:);
%             if ~(isempty(CIobj.VolDiff) && isempty(CIobj.VolCum))
%                 if length(CIobj.VolDiff)~=length(CIobj.DoseBins_org) && length(CIobj.VolCum)~=length(CIobj.DoseBins_org)
%                     warning('the (to be assigned) DoseBins_org has different number of elements from both VolDiff and VolCum, in classComplicationIndividual');
%                 end
%             end
        end
        function CIobj = set.VolDiff(CIobj,voldiff)
            if any(isnan(voldiff))
                error('there is nan element when atempting to assigning member VolDiff, in classComplicationIndividual');
            end
            CIobj.VolDiff = voldiff(:);
            if voldiff(end) ~= 0 % the last bin should be with zero volume
                CIobj.VolDiff(end+1) = 0;
            end
            if length(CIobj.VolDiff)~=length(CIobj.DoseBins_org)
                d = length(CIobj.DoseBins_org) - length(CIobj.VolDiff);
                if d>0
                    CIobj.VolDiff(end+1:end+d) =0;
                    warning('DoseBins_org has more number of elements than the (to be assigned) VolDiff, in classComplicationIndividual. VolDiff is complemented by zeros');
                else
                    warning('DoseBins_org has less number of elements than the (to be assigned) VolDiff, in classComplicationIndividual, which is complemented by expolation');
                    if length(CIobj.DoseBins_org)<2
                        dstep = 1;
                    else
                        dstep = CIobj.DoseBins_org(end)-CIobj.DoseBins_org(end-1);
                    end
                    CIobj.DoseBins_org(end:end-d) = CIobj.DoseBins_org(end):dstep:(CIobj.DoseBins_org(end)-d*dstep);
                end
            end
        end
        function CIobj = set.VolCum(CIobj,volcum)
            if any(isnan(volcum))
                error('there is nan element when atempting to assigning member VolCum, in classComplicationIndividual');
            end
            CIobj.VolCum = volcum(:);
            if volcum(end) ~= 0 % the last bin should be with zero volume
                CIobj.VolCum(end+1) = 0;
            end
            if length(CIobj.VolCum)~=length(CIobj.DoseBins_org)
                d = length(CIobj.DoseBins_org) - length(CIobj.VolCum);
                if d>0
                    warning('DoseBins_org has more number of elements than the (to be assigned) VolCum, in classComplicationIndividual, VolCum is complemented by zeros');
                    CIobj.VolCum(end+1:end+d) =0;
                else
                    warning('DoseBins_org has less number of elements than the (to be assigned) VolCum, in classComplicationIndividual, which is complemented by expolation');
                    if length(CIobj.DoseBins_org)<2
                        dstep = 1;
                    else
                        dstep = CIobj.DoseBins_org(end)-CIobj.DoseBins_org(end-1);
                    end
                    CIobj.DoseBins_org(end:end-d) = CIobj.DoseBins_org(end):dstep:(CIobj.DoseBins_org(end)-d*dstep);
                end
            end
        end
        
        function CIobj = LinearQuartraticCorrection(CIobj)
            if isempty(CIobj.Beta2Alpha)
                error('ratio of Beta to Alpha is not determined');
            end
            if CIobj.Beta2Alpha > 0
                CIobj.DoseBins_LQ = CIobj.DoseBins_org .* ( 1 + CIobj.Beta2Alpha / CIobj.FxNum * CIobj.DoseBins_org );
            else
                CIobj.DoseBins_LQ = CIobj.DoseBins_org;
            end
        end
        function CIobj = Cum2Diff(CIobj) % compute differential DVH from cumulative DVH
            CIobj.VolDiff=zeros(size(CIobj.VolCum));
            CIobj.VolDiff(1:end-1) = diff( CIobj.VolCum);
            if any(CIobj.VolDiff(:)>0) % check if it is a monotonic decreasing cumulative VolCum
                error('The cumulative VolCum is not monotonic decreasing');
            end
            CIobj.VolDiff = abs( CIobj.VolDiff);
        end
        function CIobj = Diff2Cum(CIobj)
            CIobj.VolCum = flipud( CIobj.VolDiff );
            CIobj.VolCum = cumsum( CIobj.VolCum );
            CIobj.VolCum = flipud( CIobj.VolCum );
        end
        function interpDose = DoseAtVol(CIobj,vol)
            if vol < 0
                error('Dose to be interpolated is at negative volume');
            end
            
            f = find( CIobj.VolCum >= vol );
            if isempty(f)
                interpDose = 0;
            else
                f = f(end);
                if CIobj.VolCum(f) == vol
                    f = find( CIobj.VolCum == vol );
                    interpDose = CIobj.DoseBins_LQ(f(1));
                else
                    interpDose = interp1( [CIobj.VolCum(f); CIobj.VolCum(f+1)], [CIobj.DoseBins_LQ(f); CIobj.DoseBins_LQ(f+1)], vol );
                end
            end
        end
        function interpVol = VolAtDose(CIobj,dose)
            if dose < 0
                error('Volume to be interpolated is at negative dose');
            end
            
            f = find( CIobj.DoseBins_LQ <= dose ); f=f(end); % f won't be empty by definition
            if CIobj.DoseBins_LQ(f) < CIobj.DoseBins_LQ(end) % not the last element, interpolate to get the best estimation of volume, for the last element, it is zero
                interpVol = interp1( [CIobj.DoseBins_LQ(f); CIobj.DoseBins_LQ(f+1)], [CIobj.VolCum(f); CIobj.VolCum(f+1)], dose );
            else
                interpVol = 0;
            end
        end
        function CIobj = lgn2n(CIobj)
            CIobj.nLyman = 10.^CIobj.lgn;
        end
        function CIobj = CalculateEUD(CIobj)
            % prepare differential DVHs
            if isempty(CIobj.VolDiff) % no differential data
                CIobj = CIobj.Cum2Diff();
            end
            % prepare LQ corrected doses
            if isempty(CIobj.DoseBins_LQ) % not corrected yet
                CIobj = CIobj.LinearQuartraticCorrection();
            end
            % normalize differential volume
            vol = CIobj.VolDiff;
            v = sum(vol);
            if v~=1
                vol = vol/v;
            end
            % shift the dose from left end of bin to the center
            dosebins = ( CIobj.DoseBins_LQ(1:end-1) + CIobj.DoseBins_LQ(2:end) )/2;
            % compute EUD
            [nn,dosebins] = ndgrid( CIobj.nLyman, dosebins );
            [~,vol] = ndgrid( CIobj.nLyman, vol(1:end-1) );
            dosebins = dosebins .^ (1./nn); % (di)^(1/n)
            dosebins = dosebins .* vol; % (di)^(1/n) * vi
            CIobj.EUD = (sum( dosebins, 2 )).^CIobj.nLyman;
        end
    end
    
    methods % plots
        function DiffDVHCurve_fig(CIobj,lw)
            if ~exist('ln','var')
                lw = 1;
            end
            plot(CIobj.DoseBins_org, CIobj.VolDiff,'LineWidth',lw);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('Dose'); ylabel('Volume Fraction');
        end
        function CumDVHCurve_fig(CIobj,lw)
            if ~exist('ln','var')
                lw = 1;
            end
            plot(CIobj.DoseBins_org, CIobj.VolCum,'LineWidth',lw);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('Dose'); ylabel('Volume');
        end
        function EUDCurve_n_fig(CIobj,lw)
            if ~exist('ln','var')
                lw = 1;
            end
            plot(CIobj.EUD, CIobj.lgn,'o-','LineWidth',lw);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('EUD'); ylabel('log_1_0(n)');
        end
        function EUDCurve_a_fig(CIobj,strMarker,lw)
            if ~exist('strMarker','var')
                strMarker = 'b';
            end
            if ~exist('ln','var')
                lw = 1;
            end
            
            a1 = gca;
            a2 = copyobj(a1,gcf);
            set(a2,'Color','none');
            set(a2,'Xtick',[]);
            
            plot(a1, CIobj.EUD, CIobj.lgn,strMarker,'LineWidth',lw);
            set(a2,'YAxisLocation','right');
            set(a2,'YTickLabel',num2str(CIobj.lgn(1:2:end)));
            
            set(a1,'YTickLabel',num2str(CIobj.lgn(end:-2:1)));
            set(a1,'xminortick','on','yminortick','on');
            set(a1,'box','on');
            xlabel(a1,'EUD');
            ylabel(a1,'log_1_0(a)');
%             ylabel(a2,'log_1_0(n)');
        end
    end
end