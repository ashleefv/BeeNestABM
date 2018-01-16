function bump_Indicator = bump(BeeBodyThreshold,current_Distance)
% Determine if this bee bumps one or more other bees based on the current
% positions
    bump_Indicator = zeros(size(current_Distance));
    bump_Indicator(current_Distance<=BeeBodyThreshold) = 1;
    %
    % Bees that are not in the chamber will have NaN as their
    % current_Distance and will give 0's for bumps with themselves
    %
    % Need to exclude the case when bees are finite and thus technically
    % bump themselves along the diagonal because that current_Distance == 0
    bump_Indicator = bump_Indicator - diag(diag(bump_Indicator));
    %
    bump_Indicator = sum(bump_Indicator,2);
    bump_Indicator(bump_Indicator>=1) = 1;
    %
end

