function ea_postprocess_multisource_pam(options,settings,side)

% Merge multisource PAMs using logical disjunction.
% IMPORTANT:fiberActivations here have the same pre-filtering with
% Phi_vector_max. That is why they contain the same fibers and can be
% merged directly
% By Butenko, konstantinmgtu@gmail.com

arguments
    options             % Lead-DBS options for electrode reconstruction and stimulation
    settings            % parameters for OSS-DBS simulation
    side                % 1 - rh, 2 - lh
end



% Files for the same pathway differ only by a trailing _<source> index.

% get all fiberActivations (make sure you don't have merged once before running it!)
fiberActivations = dir([settings.connectomeActivations,filesep,'sub-*']);

if side == 1
    % remove left side fiberActivations
    fiberActivations = fiberActivations(~contains({fiberActivations.name}, {'-L_','-L.mat'}));
else
    fiberActivations = fiberActivations(~contains({fiberActivations.name}, {'-R_','-R.mat'}));
end

already_merged_index = [];
for pam_idx = 1:size(fiberActivations,1)

    thisName = fiberActivations(pam_idx).name;
    % Strip a trailing _<source> before .mat. No-op when the file carries no
    % source suffix, which is what cluster postprocessing writes. The previous
    % fixed-offset chop (end-6 / end-5) ate two characters of the tract name.
    baseName = regexprep(thisName, '_(\d+)\.mat$', '.mat');
    groupKey = regexprep(thisName, '_(\d+)\.mat$', '');
    hasSuffix = ~strcmp(baseName, thisName);

    if ~ismember(pam_idx,already_merged_index)

        ftr = load([fiberActivations(pam_idx).folder,filesep,fiberActivations(pam_idx).name]);
        %[fibers_unique,ia,ic] = unique(ftr.fibers(:,4));

        fiberActivationsMerged = [fiberActivations(pam_idx).folder,filesep,baseName];

        % The MNI twin is a different fiber set, not a warped copy: it carries
        % the pristine template tract, while the native file holds only the
        % subject's prefiltered subset. Row counts therefore differ for most
        % tracts, so each space is merged against its own sources rather than
        % having a status column transplanted across spaces.
        mergeMNI = options.native ...
            && exist([settings.connectomeActivationsMNI,filesep,thisName],'file');
        if mergeMNI
            ftrMNI = load([settings.connectomeActivationsMNI,filesep,thisName]);
        end

        % search the rest to find the same pathway
        for pam_idx_2 = pam_idx:size(fiberActivations,1)
            % check if the same name except the source index
            if strcmp(groupKey, regexprep(fiberActivations(pam_idx_2).name, '_(\d+)\.mat$', ''))

                % merge here
                ftr2 = load([fiberActivations(pam_idx_2).folder,filesep,fiberActivations(pam_idx_2).name]);
                % we only need to flip 0 to 1 if necessary. Other statutes
                % must stay the same per definition!

                ftr.fibers(ftr2.fibers(:,5) == 1,5) = 1;

                if mergeMNI
                    mniSource = [settings.connectomeActivationsMNI,filesep,fiberActivations(pam_idx_2).name];
                    if exist(mniSource,'file')
                        ftr2MNI = load(mniSource);
                        ftrMNI.fibers(ftr2MNI.fibers(:,5) == 1,5) = 1;
                    end
                end

                already_merged_index = [already_merged_index, pam_idx_2];

            end
        end
        % save the merged fiberActivation
        save(fiberActivationsMerged, '-struct', 'ftr');

        if mergeMNI
            fiberActivationsMergedMNI = [settings.connectomeActivationsMNI,filesep,baseName];
            save(fiberActivationsMergedMNI, '-struct', 'ftrMNI');
        end
    end

    % remove the per-source result to avoid wrong FF imports. An unsuffixed
    % file is the merged output itself, never a source result, so it survives.
    if hasSuffix
        ea_delete([fiberActivations(pam_idx).folder,filesep,thisName])

        if options.native
            if exist([settings.connectomeActivationsMNI,filesep,thisName], 'file')
                ea_delete([settings.connectomeActivationsMNI,filesep,thisName])
            end
        end
    end

end

