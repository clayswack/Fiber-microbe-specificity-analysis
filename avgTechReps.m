function [rbFinal, rbInitial] = avgTechReps(rbFinalDouble, rbInitialDouble, initialMeta, finalMeta, taxaNames, unqFibers, nFibers, unqDonors, nDonors)
%% Objective: Average out technical replicates (from the same donor on the same fiber)
% from the relative abundance tables (initial and final timepoints)

    rbInitNoTechRepsMeta = table;
    rbFinalNoTechRepsMeta = table;
    
    
    % First average out tech reps in the final samples
    j = 1; % j indicates rows in our out data that has the technical replicates averaged out
    for i = 1:nFibers
        fiberNow = unqFibers(i);
    
        for k = 1:nDonors
            donorNow = unqDonors(k);
    
            % RBfinal
            rbFinalNow = rbFinalDouble(finalMeta.donor == donorNow & finalMeta.fiber_type == fiberNow, :);
            rbFinalNoTechReps(j,:) = mean(rbFinalNow, 'omitnan');
    
            % Save the metadata
            rbFinalNoTechRepsMeta_donors(j) = donorNow;
            rbFinalNoTechRepsMeta_fiber_type(j) = fiberNow;
            
            % Increment counter
            j = j + 1;
        end
    
    end
    
    % Assemble outputs for final samples
    rbFinalData = array2table(rbFinalNoTechReps, "VariableNames", taxaNames);
    rbFinal = table;
    rbFinal.donor = rbFinalNoTechRepsMeta_donors';
    rbFinal.fiber_type = rbFinalNoTechRepsMeta_fiber_type';
    rbFinal.time = categorical(repmat({'after'},height(rbFinalData),1));
    rbFinal = horzcat(rbFinal, rbFinalData);
    
    
    % Now average out reps in the initial samples
    j = 1;
    for k = 1:nDonors
        donorNow = unqDonors(k);
        
        % RBinitial
        rbInitNow = rbInitialDouble(initialMeta.donor == donorNow,:);
        rbInitialNoTechReps(j,:) = mean(rbInitNow, 'omitnan');
        
        % Save the metadata
        rbInitialNoTechRepsMeta_donor(j) = donorNow;
    
        % Increment counter
        j = j + 1;
    end
    
    % Assemble outputs for initial samples
    rbInitialData = array2table(rbInitialNoTechReps, "VariableNames", taxaNames);
    rbInitial = table;
    rbInitial.time = categorical(repmat({'before'},height(rbInitialNoTechReps),1));
    rbInitial.donor = rbInitialNoTechRepsMeta_donor';
    rbInitial = horzcat(rbInitial, rbInitialData);

end
