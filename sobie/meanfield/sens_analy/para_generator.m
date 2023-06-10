function [] = para_generator(name, idx)
    if idx == 1
        %% Generator k for each parameter based on 0.1-insensitive range
        % alphalst1 = linspace(0.1,4.2,83); alphalst2 = linspace(4.35,4.50,4); alphalst3 = linspace(4.60,4.60,1);
        alphalst = [1];
        % dryrlst1 = linspace(0.80,1.35,12); dryrlst2 = linspace(1.45,1.45,1);
        dryrlst = linspace(1.00,1.10,3);
        nryrlst = linspace(0.95,1.00,2);
        bcsqlst = linspace(0.90,1.00,3);
        effluxlst = linspace(1.00,1.00,1);
        refilllst = linspace(0.95,1.15,3);
        vjsrlst = linspace(0.95,1.05,3);
        vsslst = linspace(1.00,1.10,3);
        bsllst = linspace(0.95,1.05,3);
        kcsqlst = linspace(0.95,1.00,2);
        kmaxlst = [1];
        
        % the order of parameters: 
        % alpha, dryr, nryr, csq, efflux, refill, vjsr, vss, bsl, kcsq, kmax
        all_parameters = zeros(1,11);
        for i1 = 1:length(alphalst)
        for i2 = 1:length(dryrlst)
        disp(i2);
        for i3 = 1:length(nryrlst)
        for i4 = 1:length(bcsqlst)
        for i5 = 1:length(effluxlst)
        for i6 = 1:length(refilllst)
        for i7 = 1:length(vjsrlst)
        for i8 = 1:length(vsslst)
        for i9 = 1:length(bsllst)
        for i10 = 1:length(kcsqlst)
        for i11 = 1:length(kmaxlst)
        all_parameters(end+1,:) = [alphalst(i1),dryrlst(i2),nryrlst(i3),bcsqlst(i4),effluxlst(i5),refilllst(i6),...
            vjsrlst(i7),vsslst(i8),bsllst(i9),kcsqlst(i10),kmaxlst(i11)];
        end
        end
        end
        end
        end
        end
        end
        end
        end
        end
        end
        
        save(name,'all_parameters');
    end
end