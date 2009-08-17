function gpl()
%GPL  Location of full copyright information for the Spot Toolbox.

%   Copyright 2009 Ewout van den Berg and Michael P. Friedlander.
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

    parts = regexp(mfilename('fullpath'),filesep,'split');
    spotpath = [filesep fullfile(parts{1:end-2})];
    
    fprintf('\n');
    fprintf('  Full copyright information for the Spot Toolbox\n');
    fprintf('  can be found in\n\n');
    fprintf('  %s%s%s\n',spotpath,filesep,'COPYING.txt');
    fprintf('\n');
    
end % function gpl