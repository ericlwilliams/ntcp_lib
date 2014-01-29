classdef classOutcomeIndividual
    properties
        mID % patient id (MRN)
        mNameLast
        mNameMiddle
        mNameFirst
        mGender % male -- 1
        mAgeAtTx % in years
        mBMI % if no bmi, -Inf
        mKPS 
        mDistanceToChestWall 
        mDistanceToSpine
        mDosePerFx
        mDateBirth % used to compute the patient age at some time point
        mDateDeath % used to compute the survival time
        mMaxVol % volume of organ
        
        mStagePreTx % mStagePreTx before the treatment
        mDateBaseline % the date of the start point (in number)
        mDateStartTx % the date of the start point (in number)
        mDateEndTx % the date of end of Tx
        mDateLastFollowUp % date of last evaluation
        mDateComp % the date of complication occurance (in number)
        mDateComp2 % the date of complication occurance (in number)
        mCompGrade % patient complication grade (w/in 6 months)
        mMaxToxGrade % maximum toxicity grade
        mDateLastFollowup % the date of the last followup (in number)
        mDateRelapse % the date of relapse (in number)
        mFlgCensor % =1, censored (e.g. >= 3Grade RP, or death)
        mFlgCensor2 % =1, censored (e.g. >= 3Grade RP, or death)
        mDateCensor % date of censor
        mFlgNFz % 1- in Nfz
        
        mDoseRx % dose prescription
        mDoseTx % dose delivered (treated)
        mDoseBins_org
        mVolDiff
        mVolCum
        mBeta2Alpha
        % USC parameters
        mDoseTrans
        mAlphaUSC
        mDzeroUSC
        mDqUSC
        
        
        mBeta2AlphaCorrection
        mDoseBins_LQ
        mFxNum

        mHeartTargetOL
        
        mLymanN % parameter n in Lyman model
        mEUD % equivalent uniform dose
    end
    
    methods % operations
        function EPIobj = classOutcomeIndividual()
        end
        function EPIobj = set.mID(EPIobj,ptid)
            EPIobj.mID = num2str(ptid);
        end
        function EPIobj = set.mLymanN(EPIobj,n)
            EPIobj.mLymanN = n(:);
        end
        function EPIobj = set.mDoseBins_org(EPIobj,dosebins)
            if any(isnan(dosebins))
                error('there is nan element when atempting to assigning member mDoseBins_org, in classComplicationIndividual');
            end
            EPIobj.mDoseBins_org = dosebins(:);
%             if ~(isempty(EPIobj.mVolDiff) && isempty(EPIobj.mVolCum))
%                 if length(EPIobj.mVolDiff)~=length(EPIobj.mDoseBins_org) && length(EPIobj.mVolCum)~=length(EPIobj.mDoseBins_org)
%                     warning('the (to be assigned) mDoseBins_org has different number of elements from both mVolDiff and mVolCum, in classComplicationIndividual');
%                 end
%             end
        end
        function EPIobj = set.mVolDiff(EPIobj,mVolDiff)
            if any(isnan(mVolDiff))
                error('there is nan element when atempting to assigning member mVolDiff, in classComplicationIndividual');
            end
            EPIobj.mVolDiff = mVolDiff(:);
            if mVolDiff(end) ~= 0 % the last bin should be with zero volume
                EPIobj.mVolDiff(end+1) = 0;
            end
            if length(EPIobj.mVolDiff)~=length(EPIobj.mDoseBins_org)
                d = length(EPIobj.mDoseBins_org) - length(EPIobj.mVolDiff);
                if d>0
                    EPIobj.mVolDiff(end+1:end+d) =0;
                    warning('mDoseBins_org has more number of elements than the (to be assigned) mVolDiff, in classComplicationIndividual. mVolDiff is complemented by zeros');
                else
                    warning('mDoseBins_org has less number of elements than the (to be assigned) mVolDiff, in classComplicationIndividual, which is complemented by expolation');
                    if length(EPIobj.mDoseBins_org)<2
                        dstep = 1;
                    else
                        dstep = EPIobj.mDoseBins_org(end)-EPIobj.mDoseBins_org(end-1);
                    end
                    EPIobj.mDoseBins_org(end:end-d) = EPIobj.mDoseBins_org(end):dstep:(EPIobj.mDoseBins_org(end)-d*dstep);
                end
            end
        end
        function EPIobj = set.mVolCum(EPIobj,mVolCum)
            if any(isnan(mVolCum))
                error('there is nan element when atempting to assigning member mVolCum, in classComplicationIndividual');
            end
            EPIobj.mVolCum = mVolCum(:);
            % tmp
            if mVolCum(end) ~= 0 % the last bin should be with zero volume
                EPIobj.mVolCum(end+1) = 0;
                %EPIobj.mVolCum(1) = [];
            end
            if length(EPIobj.mVolCum)~=length(EPIobj.mDoseBins_org)
                d = length(EPIobj.mDoseBins_org) - length(EPIobj.mVolCum);
                if d>0
                    warning('mDoseBins_org has more number of elements than the (to be assigned) mVolCum, in classComplicationIndividual, mVolCum is complemented by zeros');
                    EPIobj.mVolCum(end+1:end+d) =0;
                else
                    warning('mDoseBins_org has less number of elements than the (to be assigned) mVolCum, in classComplicationIndividual, which is complemented by expolation');
                    if length(EPIobj.mDoseBins_org)<2
                        dstep = 1;
                    else
                        dstep = EPIobj.mDoseBins_org(end)-EPIobj.mDoseBins_org(end-1);
                    end
                    EPIobj.mDoseBins_org(end:end-d) = EPIobj.mDoseBins_org(end):dstep:(EPIobj.mDoseBins_org(end)-d*dstep);
                end
            end
        end
        function EPIobj = set.mBeta2Alpha(EPIobj,beta2alpha)
            if length(beta2alpha)>1
                warning('There are more than one Beta/Alpha available, use the first one');
            end
            EPIobj.mBeta2Alpha = beta2alpha(1);
            EPIobj = EPIobj.fLinearQuartraticCorrection();
        end
        function EPIobj = fUndoLinearQuadratic(EPIobj,alpha2beta)
            % solve quadratic equation
            quad_a = (1/EPIobj.mFxNum);
            quad_b = alpha2beta;
            quad_c = -(alpha2beta + 2).*EPIobj.mDoseBins_org;
            
            %quad_sols = zeros(2,1);
            quad_d = sqrt(quad_b^2 - 4*quad_a*quad_c);
            quad_sols1 = (-quad_b+quad_d)/(2*quad_a);
            quad_sols2 = (-quad_b-quad_d)/(2*quad_a);
            
            % first solution always positive, 
            if sum(quad_sols1<0) >0 || sum(quad_sols2<0)==0,
                error(['Problem with quadratic solution in fUndoLinearQuadratic']);
            end

            EPIobj.mDoseBins_org = quad_sols1;
            
        end
            
        function EPIobj = fLinearQuartraticCorrection(EPIobj)
            if isempty(EPIobj.mBeta2Alpha)
                error('ratio of Beta to Alpha is not determined');
            end
            if ((EPIobj.mBeta2Alpha > 0) && ~isnan(EPIobj.mFxNum) && ~isempty(EPIobj.mBeta2AlphaCorrection))
                if (strcmp(EPIobj.mBeta2AlphaCorrection,'BED'))
                    EPIobj.mDoseBins_LQ = EPIobj.mDoseBins_org .* ( 1 + ((EPIobj.mDoseBins_org/EPIobj.mFxNum)*EPIobj.mBeta2Alpha));
                elseif (strcmp(EPIobj.mBeta2AlphaCorrection,'NTD'))
                    EPIobj.mDoseBins_LQ = EPIobj.mDoseBins_org .* ...
                        (((1/EPIobj.mBeta2Alpha) + (EPIobj.mDoseBins_org/EPIobj.mFxNum))/...
                        ((1/EPIobj.mBeta2Alpha) + 2));
                elseif (strcmp(EPIobj.mBeta2AlphaCorrection,'USCBED'))
                    if isempty(EPIobj.mDoseTrans)
                        error('D_{trans} for USCBED not calculated');
                    end
                    use_lq = logical([(EPIobj.mDoseBins_org/EPIobj.mFxNum) <= EPIobj.mDoseTrans]);
                    EPIobj.mDoseBins_LQ = inf(size(EPIobj.mDoseBins_org));
                    
                    %lq
                   EPIobj.mDoseBins_LQ(use_lq) =  EPIobj.mDoseBins_org(use_lq) .* ( 1 + ((EPIobj.mDoseBins_org(use_lq)/EPIobj.mFxNum)*EPIobj.mBeta2Alpha));
                    
                   %usc bed
                    EPIobj.mDoseBins_LQ(~use_lq) = (1/(EPIobj.mAlphaUSC*EPIobj.mDzeroUSC))*...
                            (EPIobj.mDoseBins_org(~use_lq) - (EPIobj.mFxNum*EPIobj.mDqUSC));
                                       
                end
            else % if mFxNum is empty, don't do correction
                EPIobj.mDoseBins_LQ = EPIobj.mDoseBins_org;
            end
        end
        
        function [EPIobj,num_lq,num_bins,full_lq] = fUSCCorrection(EPIobj,dq,alfad0,alfa2beta)
                    
            dt = (2*dq)/(1-(alfad0));

            use_lq = logical([(EPIobj.mDoseBins_org/EPIobj.mFxNum) <= dt]);
            
            num_lq = sum(use_lq);
            num_bins= length(use_lq);

            if isequal(num_lq,num_bins)
                full_lq = 1;
            else
                full_lq = 0;
            end
                
            EPIobj.mDoseBins_LQ = inf(size(EPIobj.mDoseBins_org));
            
            %lq
                EPIobj.mDoseBins_LQ(use_lq) =  EPIobj.mDoseBins_org(use_lq) .* ( 1 + ((EPIobj.mDoseBins_org(use_lq)/EPIobj.mFxNum)*(1/alfa2beta)));
            %usc bed
            EPIobj.mDoseBins_LQ(~use_lq) = (1/(alfad0))*...
                (EPIobj.mDoseBins_org(~use_lq) - (EPIobj.mFxNum*dq));
            
            
        end
        
        
        
        function EPIobj = fCum2Diff(EPIobj) % compute differential DVH from cumulative DVH
            EPIobj.mVolDiff=zeros(size(EPIobj.mVolCum));
            EPIobj.mVolDiff(1:end-1) = diff( EPIobj.mVolCum);
            if any(EPIobj.mVolDiff(:)>0) % check if it is a monotonic decreasing cumulative mVolCum
                error('The cumulative mVolCum is not monotonic decreasing');
            end
            EPIobj.mVolDiff = abs( EPIobj.mVolDiff);
        end
        function EPIobj = fDiff2Cum(EPIobj)
            %EPIobj.mVolCum = flipud( EPIobj.mVolDiff );
            tmp_vol = EPIobj.mVolDiff;
            tmp_vol = flipud(EPIobj.mVolDiff);
            tmp_vol = cumsum(tmp_vol);
            EPIobj.mVolCum = flipud( tmp_vol);
            %EPIobj.mVolCum = flipud( EPIobj.mVolDiff);
            %EPIobj.mVolCum = cumsum( EPIobj.mVolCum );
            %EPIobj.mVolCum = flipud( EPIobj.mVolCum );
        end
        function interpDose = fDoseAtVol(EPIobj,vol)
            if vol < 0
                error('Dose to be interpolated is at negative volume');
            end
            
            f = find( EPIobj.mVolCum >= vol );
            if isempty(f)
                interpDose = 0;
            else
                f = f(end);
                if EPIobj.mVolCum(f) == vol
                    f = find( EPIobj.mVolCum == vol );
                    interpDose = EPIobj.mDoseBins_LQ(f(1));
                else
                    interpDose = interp1( [EPIobj.mVolCum(f); EPIobj.mVolCum(f+1)], [EPIobj.mDoseBins_LQ(f); EPIobj.mDoseBins_LQ(f+1)], vol );
                end
            end
        end
        function interpVol = fVolAtDose(EPIobj,dose)
            if dose < 0
                error('Volume to be interpolated is at negative dose');
            end
            
            f = find( EPIobj.mDoseBins_LQ <= dose ); f=f(end); % f won't be empty by definition
            if EPIobj.mDoseBins_LQ(f) < EPIobj.mDoseBins_LQ(end) % not the last element, interpolate to get the best estimation of volume, for the last element, it is zero
                interpVol = interp1( [EPIobj.mDoseBins_LQ(f); EPIobj.mDoseBins_LQ(f+1)], [EPIobj.mVolCum(f); EPIobj.mVolCum(f+1)], dose );
            else
                interpVol = 0;
            end
        end
      
        function EPIobj = fCalculateEUD(EPIobj)
            % prepare differential DVHs
            if isempty(EPIobj.mVolDiff) % no differential data
                EPIobj = EPIobj.fCum2Diff();
            end
            % prepare LQ corrected doses
            if isempty(EPIobj.mDoseBins_LQ) % not corrected yet
                EPIobj = EPIobj.fLinearQuartraticCorrection();
            end
            % normalize differential volume
            vol = EPIobj.mVolDiff;
            v = sum(vol);
            if v~=1
                vol = vol/v;
            end
            % shift the dose from left end of bin to the center
            dosebins = ( EPIobj.mDoseBins_LQ(1:end-1) + EPIobj.mDoseBins_LQ(2:end) )/2;
            % compute EUD
            [nn,dosebins] = ndgrid( EPIobj.mLymanN, dosebins );
            [~,vol] = ndgrid( EPIobj.mLymanN, vol(1:end-1) );
            dosebins = dosebins .^ (1./nn); % (di)^(1/n)
            dosebins = dosebins .* vol; % (di)^(1/n) * vi
            EPIobj.mEUD = (sum( dosebins, 2 )).^EPIobj.mLymanN;
        end
     end
    
    methods % plots
        function fDiffDVHCurve_fig(EPIobj,lw)
            if ~exist('ln','var')
                lw = 1;
            end
            plot(EPIobj.mDoseBins_org, EPIobj.mVolDiff,'LineWidth',lw);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('Dose'); ylabel('Volume Fraction');
        end
        function fCumDVHCurve_fig(EPIobj,lw)
            if ~exist('ln','var')
                lw = 1;
            end
            plot(EPIobj.mDoseBins_org, EPIobj.mVolCum,'LineWidth',lw);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('Dose'); ylabel('Volume');
        end
        function fEUDCurve_log10n_fig(EPIobj,lw)
            if ~exist('ln','var')
                lw = 1;
            end
            plot(EPIobj.mEUD, log10(EPIobj.mLymanN),'o-','LineWidth',lw);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('EUD'); ylabel('log_1_0(n)');
        end
        function fEUDCurve_log10a_fig(EPIobj,strMarker,lw)
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
            
            plot(a1, EPIobj.mEUD, log10(EPIobj.mLymanN),strMarker,'LineWidth',lw);
            set(a2,'YAxisLocation','right');
            set(a2,'YTickLabel',num2str(log10(EPIobj.mLymanN(1:2:end))));
            
            set(a1,'YTickLabel',num2str(log10(EPIobj.mLymanN(end:-2:1))));
            set(a1,'xminortick','on','yminortick','on');
            set(a1,'box','on');
            xlabel(a1,'mEUD');
            ylabel(a1,'log_1_0(a)');
%             ylabel(a2,'log_1_0(n)');
        end
    end
end