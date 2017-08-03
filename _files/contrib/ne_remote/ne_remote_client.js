/*

    NeuroElf remote web interface
    JavaScript component

% Version:  v0.9c
% Build:    12012418
% Date:     Jan-24 2012, 6:07 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2012, Jochen Weber
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

*/

// global document variable
var NE_REMOTE = this;

// background requests
var bgrmode = true;

// controls (so we don't need to get them over and over again...)
var xcrd = ycrd = zcrd = null;

// current position
var cpos = new Array(0, 0, 0);

// last state
var laststate = '';

// tc-plot data
var tcplotx, tcplotdata, tcplotsyms;
var tcplotsvg = d3.select("#tcplotsvg");

function NE_ObjType(o) {

    if (o===null) return '[object null]';
    return Object.prototype.toString.call(o);
}

function NE_Sleep(delay) {
    var start = new Date().getTime();
    while (new Date().getTime() < start + delay);
}

function NE_HTTPRequest() {

    try {
        return new XMLHttpRequest();
    } catch (error) { }

    try {
        return new ActiveXObject("Msxml2.XMLHTTP");
    } catch (error) { }

    try {
        return new ActiveXObject("Microsoft.XMLHTTP");
    } catch (error) { }

    throw new Error("Could not create HTTP request object.");
}

function NE_HTTPGet(getURI, listener) {

    var request = NE_HTTPRequest();

    if (arguments.length > 1) {
        request.onreadystatechange = function () {
                if (request.readyState == 4) {
                    if (request.responseText.match(/^session_timeout/)) {
                        var ddate = new Date();
                        document.location.href = 'timeout.html';
                    } else
                        listener(request.responseText);
                }
            }
        request.open("GET", getURI, true);
    } else
        request.open("GET", getURI, false);
    request.send(null);

    if (arguments.length < 2) {
        var rsc = 0;
        while (request.readyState < 4) {

            rsc++;
            if (rsc > 1500) { break; }
            NE_sleep(10);
        }

        if (request.readyState < 4 || request.status != 200) {
            throw new Error("Request timed out or returned with an error.");
        }

        return request.responseText;
    } else {
        return request;
    }
}

function NE_OpenFile(file) {

    NE_RemoteCMD("openfile%25" + file, function(rstate) {
        if (rstate.match(/OK/)) {
            NE_ReadUIState();
            NE_UpdateSliceImgs();
        }
    });
}

function NE_ReadDir(folder) {

    if (arguments.length < 1 || !folder) {
        var folder = '/';
        var u = document.URL;
        var ua = u.split('/');
        u = ua[ua.length-1];
        ua = u.split('?');
        u = ua[1];
        ua = u.split('&');
        u = '';
        for (var uc = 0; uc < ua.length; uc++) {
            if (ua[uc].match(/^dir/)) {
                u = ua[uc];
                ua = u.split('=');
                break;
            }
        }
        if ((u != '') && (ua.length > 0))
            folder = ua[1];
	}

	NE_RemoteCMD("dir%25" + folder, function(rstate) {
        var listc = rstate.split('\\n');
        var liste = document.getElementById('LB_NeuroElf_DirList');
        var ih = '<li><b>Folder contents of: ' +
			listc[0].substring(0, listc[0].length - 1) + '</b></li>';
        for (var lc = 1; lc < listc.length; lc++) {
            if (listc[lc].match(/\/$/))
                ih = ih + '<li><a href="#" onclick="NE_ReadDir(\'' +
                    listc[0] + listc[lc] + '\');return false;">' +
                    listc[lc] + '</a></li>';
            else
                ih = ih + '<li><a href="#" onclick="NE_OpenFile(\'' +
                    listc[0] + listc[lc] + '\');return false;">' +
                    listc[lc] + '</a></li>';
        }
        liste.innerHTML = ih;
    });
}

function NE_RemoteCMD(cmdPart, listener) {

    if (arguments.length > 1)
        return NE_HTTPGet("?cmd=" + cmdPart, listener);
    else
        return NE_HTTPGet("?cmd=" + cmdPart);
}

function NE_ReloadImg(imgid) {

    var imgobj = document.getElementById(imgid);
    var imgsrc = imgobj.src;
    var pos = imgsrc.indexOf('&');
    if (pos >= 0) {
        imgsrc = imgsrc.substr(0, pos);
    }
    var date = new Date();
    imgobj.src = imgsrc + '&update=' + date.getTime();
    return false;
}

function NE_ReloadTCPlot(plotdata) {

    if (arguments.length < 1) {
        var plotdata = NE_RemoteCMD("tcplot", function(rplotdata) { NE_ReloadTCPlot(rplotdata); });
        return;
    }
    d3.csv.parse(plotdata, function(data) {
        tcplotsyms = d3.nest().key(function(d) { return d.symbol; }).entries(tcplotdata = data);
        tcplotsyms.forEach(function(s) {
            s.values.forEach(function(d) { d.index = +d.index; d.value = +d.value; });
            s.maxvalue = d3.max(s.values, function(d) { return d.value; });
        });
    });

    var g = tcplotsvg.data(tcplotsyms).enter().append("g").attr("class", "symbol");
}

function NE_UpdateSliceImgs() {
    NE_ReloadImg("SagSlice");
    NE_ReloadImg("CorSlice");
    NE_ReloadImg("TraSlice");
    NE_ReloadTCPlot();
}

function NE_ReadCPos(rstate, upui) {

    var cposstrA = rstate.split("%");
    cpos[0] = parseInt(cposstrA[0]);
    cpos[1] = parseInt(cposstrA[1]);
    cpos[2] = parseInt(cposstrA[2]);

    xcrd.value = cpos[0].toString();
    ycrd.value = cpos[1].toString();
    zcrd.value = cpos[2].toString();

    if (arguments.length > 1 && upui) {
        NE_UpdateSliceImgs();
    }
}

function NE_GetCPos(update, upui) {

    if (arguments.length > 0) {
        NE_RemoteCMD("cpos" + "%25" + cpos[0].toString() + "%25" + cpos[1].toString() + "%25" + cpos[2].toString(), function(rstate) { NE_ReadCPos(rstate, upui); });
    } else {
        NE_RemoteCMD("cpos", function(rstate) { NE_ReadCPos(rstate, upui); });
    }
}

function NE_AbsoluteObjectPosition(o) {

    var l = 0;
    var t = 0;
    if (o.offsetParent) {
        do {
            l += o.offsetLeft;
            t += o.offsetTop;
        } while (o = o.offsetParent);
    }
    return [l, t];
}

function NE_AbsoluteMousePosition(e) {

    var posx = 0;
    var posy = 0;
    var e = e ? e : window.event;
    if (e.pageX || e.pageY) 	{
        posx = e.pageX;
        posy = e.pageY;
    } else if (e.clientX || e.clientY) 	{
        posx = e.clientX + document.body.scrollLeft + document.documentElement.scrollLeft;
        posy = e.clientY + document.body.scrollTop + document.documentElement.scrollTop;
	}
    return [posx, posy];
}

function NE_RelativeMousePosition(e, mobj) {

    var mpos = NE_AbsoluteMousePosition(e);
    var opos = NE_AbsoluteObjectPosition(mobj);
    var rposx = mpos[0] - opos[0];
    var rposy = mpos[1] - opos[1];
    return [rposx, rposy];
}

function NE_BlurOnReturn(e, uic, upui) {

    var e = e ? e : window.event;
    var upui = upui ? upui : false;
    if (e.keyCode == 13 || e.keyCode == 10) {
        NE_OnBlur(uic, upui);
        return true;
    } else
        return false;
}

function NE_Click(uic, upui, bgclick, cargs) {

    var cargstr = '';
    if (arguments.length > 3 && cargs.length > 0) {
        cargstr = "%25remote";
        for (var cac = 0; cac < cargs.length; cac++) {
            cargstr = cargstr + "%25" + cargs[cac];
        }
    }
    if (arguments.length < 3)
        var bgclick = false;
    var bgcstring = '';
    if (bgclick)
        bgcstring = '%25background';
    if (bgrmode) {
        if ((arguments.length > 1) && upui)
            NE_RemoteCMD("click%25" + uic + bgcstring + cargstr, function(rstate) {
                NE_ReadUIState(rstate);
                NE_UpdateSliceImgs();
            });
        else
            NE_RemoteCMD("click%25" + uic, function(rstate) { return rstate; });
    } else {
        var state = NE_RemoteCMD("click%25" + uic + cargstr);
        if ((arguments.length > 1) && upui) {
            NE_ReadUIState(state);
            NE_UpdateSliceImgs();
        }
    }
}

function NE_OnBlur(uic, upui) {

    var uice = document.getElementById(uic);
    if (!uice)
        return;

    var state = '';
    var uicstring = '';
    if (uic.match(/^ED/)) {
        uicstring = uice.value;

    } else if (uic.match(/^TX/)) {
        if (uict.match(/input/i))
            uicstring = uice.value;
        else
            uicstring = uice.innerHTML.replace(/\<br\s*\/?\>/g, '\n');
    }

    uicstring = "setstring%25" + uic + "%25" + encodeURIComponent(uicstring);
    if (arguments.length > 0)
        uicstring = uicstring + "%251";

    state = NE_RemoteCMD(uicstring);

    if (arguments.length > 1 && upui)
        NE_ReadUIState(state);
}

function NE_Select(uic, upui) {

    var uice = document.getElementById(uic);
    if (!uice)
        return;

    var newidx = new Array();
    var nii = 0;
    for (var gic = 0; gic < uice.options.length; gic++) {
        if (uice.options[gic].selected) {
            newidx[nii] = gic + 1;
            nii++;
        }
    }

    var reqstr = "setvalue%25" + uic + "%25";
    if (newidx.length > 0) {
        var reqidx = newidx[0].toString();
        for (var nic = 1; nic < newidx.length; nic++) {
            reqidx = reqidx + "%2c" + newidx[nic].toString();
        }
        reqstr = reqstr + reqidx;
    } else {
        reqstr = reqstr + "0";
    }
    if (bgrmode) {
        NE_RemoteCMD(reqstr, function(rstate) {
                NE_ReadUIState(rstate);
                NE_UpdateSliceImgs();
            });
    } else {
        NE_RemoteCMD(reqstr);

        if (arguments.length > 1 && upui) {

            NE_UpdateSliceImgs();
            NE_ReadUIState();
        }
    }
}

function NE_SetValue(uic, upui) {

    var uice = document.getElementById(uic);
    if (!uice)
        return;

    var state = '';
    if (uice.checked)
        state = NE_RemoteCMD("setvalue%25" + uic + "%251");
    else
        state = NE_RemoteCMD("setvalue%25" + uic + "%250");

    if (upui) {
        NE_UpdateSliceImgs();
        NE_ReadUIState();
    }
}

function NE_SetVolume(uic) {

    var uice = document.getElementById(uic);
    if (!uice)
        return;

    var state = NE_RemoteCMD("setvolume%25" + uice.value);
}

function NE_SetVolumeOnReturn(e, uic, upui) {

    var e = e ? e : window.event;
    var upui = upui ? upui : false;
    if (e.keyCode == 13 || e.keyCode == 10) {
        NE_SetVolume(uic, upui);
        return true;
    } else
        return false;
}

function NE_UpdateSlicePos(e, imgid) {

    var rpos = NE_RelativeMousePosition(e, document.getElementById(imgid));

    switch (imgid) {

        case 'SagSlice':
            cpos[1] = 128 - rpos[0];
            cpos[2] = 128 - rpos[1];
            break;
        case 'CorSlice':
            cpos[0] = rpos[0] - 128;
            cpos[2] = 128 - rpos[1];
            break;
        case 'TraSlice':
            cpos[0] = rpos[0] - 128;
            cpos[1] = 128 - rpos[1];
            break;
    }
    NE_GetCPos(true, true);
}

function NE_SetListContent(o, list, sellist, multi) {

    o.options.length = 0;
    var lc = 0;
    for (lc = 0; lc < list.length; lc++)
        o.add(new Option(list[lc], lc));
    for (lc = 0; lc < sellist.length; lc++) {
        if (sellist[lc].match(/[0-9]+/)) {
            o.options[parseInt(sellist[lc]) - 1].selected = true;
            if (!multi)
                break;
        }
    }
}

function NE_ReadUIState(state) {

    if (arguments.length < 1) {
        var state = NE_RemoteCMD("uistate", function(rstate) { NE_ReadUIState(rstate); });
        return;
    }
    var stateA = state.split("%");

    for (var ac = 0; ac < stateA.length; ac += 4) {

        var uic = stateA[ac];
        var uice = document.getElementById(uic);
        if (uice) {

            var uict = NE_ObjType(uice);

            if (stateA[ac+3].match(/on/))
                uice.disabled = false;
            else
                uice.disabled = true;

            if (uic.match(/^CB/)) {
                if (stateA[ac+2].match(/1/))
                    uice.checked = true;
                else
                    uice.checked = false;

            } else if (uic.match(/^DD/)) {
                NE_SetListContent(uice, stateA[ac+1].split(/\\n/), stateA[ac+2].split(/\,/), false);

            } else if (uic.match(/^ED/)) {
                uice.value = stateA[ac+1].replace(/\\n/g, '\n');

            } else if (uic.match(/^LB/)) {
                NE_SetListContent(uice, stateA[ac+1].split(/\\n/), stateA[ac+2].split(/\,/), true);

            } else if (uic.match(/^TX/)) {
                if (uict.match(/input/i))
                    uice.value = stateA[ac+1].replace(/\\n/g, '\n');
                else
                    uice.innerHTML = stateA[ac+1].replace(/\\n/g, '<br />');

            }
        }
    }
}

function NE_RemoteInit() {

    xcrd = document.getElementById('ED_NeuroElf_TALX');
    ycrd = document.getElementById('ED_NeuroElf_TALY');
    zcrd = document.getElementById('ED_NeuroElf_TALZ');
    NE_ReadUIState();
    NE_ReadDir('/');
}
