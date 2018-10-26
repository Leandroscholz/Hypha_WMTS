function summary=resultsSummary(resultsCellArray,Deltax)
% gets summary results from the solution cell array
% inputs 
%       resultsCellArray, result output from mainWithEventAndLoop.m
%       Deltax  -  reactor length, dm

% get the number of runs 
nRuns = length(resultsCellArray);
maxN = length(resultsCellArray{end});

summary = zeros(nRuns,maxN);

for i=1:nRuns
    
    % gets time, nutrient and vesicle concentrations but discards tip-tank
    % reactor length value from the last element in the resultsCellArray.
    % Then stores runResults into summary
    runResults = resultsCellArray{i}(end,1:end-1);
    
    % computes the number of reactors N from the size of runResults,
    % then computes final length of the run and finally stores it to the
    % last column of the summary matrix
    nReactors = (length(runResults)-1)/2; 
    finalLength = nReactors*Deltax; 
    summary(i,end) = finalLength;
    
    % Then stores runResults into summary
    summary(i,1:length(runResults)) = runResults;
end

figure(1)
plot(summary(:,1),summary(:,end)*1e5);
ylabel('Length of the hypha, \mu m')
xlabel('Simulation time, h')

