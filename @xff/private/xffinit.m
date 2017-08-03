function xffinit(xffver)
%XFFINIT  Perform initialization of XFF class singleton/factory.
%   This is an internally called function and should not be used outside
%   of the xff class constructor (first call).
%
%   See also XFF

% Version:  v1.1
% Build:    16061411
% Date:     Jun-14 2016, 11:58 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of Columbia University nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% global singleton
global xffsngl;

% init storage
xffsngl = struct( ...
    'BFF',  [], ...
    'CONF', [], ...
    'EXT',  [], ...
    'FF',   [], ...
    'FM',   [], ...
    'INIT', false, ...
    'LAST', [], ...
    'MAG',  [], ...
    'NELF', neuroelf, ...
    'OBJS', {cell(0, 4)}, ...
    'TFF',  [], ...
    'VER',  xffver);

% on initialization, also init configuration
xffconf = struct;

% set last object
xffconf.last = [0, -1];

% loaded objects of private functions
xffconf.loadedobjs = struct;

% open sliceable files in (a valid) GUI
xffconf.loadingui = false;

% load gzip-ed files
xffconf.loadgziped = true;

% reloadsame controls whether the same file (name) should be
% opened again on xff(filename) calls or a handle to the
% already open object should be returned instead
xffconf.reloadsame = true;

% settings are loaded from $config/xff.ini and control
% some of the functions (please check the contents of this file
% every now and then)
try
    bsi = xini([neuroelf_path('config') '/xff.ini'], 'convert');
    xffconf.settings = bsi.GetComplete;
    bsi.Release;
catch xfferror
    rethrow(xfferror);
end

% override empty or non-existing GZip.TempDir
if isempty(xffconf.settings.GZip.TempDir) || ...
    exist(xffconf.settings.GZip.TempDir, 'dir') ~= 7
    xffconf.settings.GZip.TempDir = tempdir;
end

% type later contains one field (1x1 struct) per filetype
xffconf.type = struct;

% update later contains one field (1x1 logical) per filetype
% if this field is set to true and a method with filename
% OBJ_Update(obj, F, S, V) exists, the method is called after
% calls to subsasgn(obj, S, V)
xffconf.update = struct;

% get file formats
cachefile = [neuroelf_path('cache') filesep 'xffcache.mat'];
xffsngl.cachefile = cachefile;
if exist(cachefile, 'file') == 2
    try
        cachemat = load(cachefile);
        cache = cachemat.cache;
        if isfield(cache, 'CONF') && isfield(cache, 'FF') && ...
            isfield(cache, 'VER') && strcmp(cache.VER, xffver) && ...
            isstruct(cache.CONF) && isfield(cache.CONF, 'settings') && ...
            isstruct(cache.CONF.settings) && isfield(cache.CONF.settings, 'GZip') && ...
            isstruct(cache.CONF.settings.GZip) && isfield(cache.CONF.settings.GZip, 'TempDir') && ...
            exist(cache.CONF.settings.GZip.TempDir, 'dir') > 0
            xffconf = cache.CONF;
            xffsngl.FF = cache.FF;
            cachefile = '';
        end
    catch xfferror
        neuroelf_lasterr(xfferror);
    end
end
if ~isempty(cachefile)
    [xffsngl.FF, xffconf] = xffformats(neuroelf_path('formats'), xffconf);
    try
        cache = struct('CONF', xffconf, 'FF', xffsngl.FF, 'VER', xffver);
        save(cachefile, 'cache');
        xffconf = cache.CONF;
    catch xfferror
        neuroelf_lasterr(xfferror);
    end
end

% hdr.assumeflipped is used when the HDR is read to determine the
% convention state if the file itself does not contain this
xffconf.type.hdr.assumeflipped = true;

% vmr.autov16 decides whether V16 files are automatically loaded/saved
xffconf.type.vmr.autov16 = true;

% get file format specific methods
xffsngl.FM = xffmethods(xffsngl.FF.extensions);

% fill main object
xffsngl.BFF  = xffsngl.FF.bff;
xffsngl.CONF = xffconf;
xffsngl.EXT  = xffsngl.FF.extensions;
%xffsngl.FF  = already filled;
%xffsngl.FM  = already filled;
xffsngl.INIT = true; % necessary to not trap next call
xffroot = xff();
xffsngl.LAST = xffroot;
xffsngl.MAG  = xffsngl.FF.magic;
xffsngl.OBJS = {'', 'ROOT', xffroot.L, xffroot};
xffsngl.TFF  = xffsngl.FF.tff;

% override transio from settings
tiotypes = xffconf.settings.Behavior.TransIOTypes;
for orfc = 1:numel(tiotypes)
    root_TransIOSize(xffroot, tiotypes{orfc}, 4096);
end
