classdef classDVHCurve
    properties
        PatientID
        DoseBins_org
        VolDiff
        VolCum
        
        FxNum=0;
        Beta2Alpha
        DoseBins_LQ
        
    end
    methods
        function DVHobj = classDVHCurve()
        end
        function DVHobj = set.VolDiff(DVHobj,voldiff)
            DVHobj.VolDiff = voldiff(:);
        end
        function DVHobj = set.VolCum(DVHobj,volcum)
            DVHobj.VolCum = volcum(:);
        end
        function DVHobj = LinearQuartraticCorrection(DVHobj)
            if DVHobj.Beta2Alpha > 0
                DVHobj.DoseBins_LQ = DVHobj.DoseBins_org .* ( 1 + DVHobj.Beta2Alpha / DVHobj.FxNum * DVHobj.DoseBins_org );
            else
                DVHobj.DoseBins_LQ = DVHobj.DoseBins_org;
            end
        end
        function DVHobj = Cum2Diff(DVHobj) % compute differential DVH from cumulative DVH
            DVHobj.VolDiff=zeros(size(DVHobj.VolCum));
            DVHobj.VolDiff(1:end-1) = diff( DVHobj.VolCum);
            if any(DVHobj.VolDiff(:)>0) % check if it is a monotonic decreasing cumulative VolCum
                error('The cumulative VolCum is not monotonic decreasing');
            end
            DVHobj.VolDiff = abs( DVHobj.VolDiff);
        end
        function DVHobj = Diff2Cum(DVHobj)
            DVHobj.VolCum = flipud( DVHobj.VolDiff );
            DVHobj.VolCum = cumsum( DVHobj.VolCum );
            DVHobj.VolCum = flipud( DVHobj.VolCum );
        end
        function interpDose = DoseAtVol(DVHobj,vol)
            if vol < 0
                error('Dose to be interpolated is at negative volume');
            end
            
            f = find( DVHobj.VolCum >= vol );
            if isempty(f)
                interpDose = 0;
            else
                f = f(end);
                if DVHobj.VolCum(f) == vol
                    f = find( DVHobj.VolCum == vol );
                    interpDose = DVHobj.DoseBins_LQ(f(1));
                else
                    interpDose = interp1( [DVHobj.VolCum(f); DVHobj.VolCum(f+1)], [DVHobj.DoseBins_LQ(f); DVHobj.DoseBins_LQ(f+1)], vol );
                end
            end
        end
        function interpVol = VolAtDose(DVHobj,dose)
            if dose < 0
                error('Volume to be interpolated is at negative dose');
            end
            
            f = find( DVHobj.DoseBins_LQ <= dose ); f=f(end); % f won't be empty by definition
            if DVHobj.DoseBins_LQ(f) < DVHobj.DoseBins_LQ(end) % not the last element, interpolate to get the best estimation of volume, for the last element, it is zero
                interpVol = interp1( [DVHobj.DoseBins_LQ(f); DVHobj.DoseBins_LQ(f+1)], [DVHobj.VolCum(f); DVHobj.VolCum(f+1)], dose );
            else
                interpVol = 0;
            end
        end
    end
end