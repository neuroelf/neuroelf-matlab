// NeuroElf - remote server
//
// testing only so far!
//
// version v0.9d, snap 14081718
//
// example URLs:
//
// http://SERVER:PORT/remote?cmd=MATLAB_COMMAND
//
// simply passes on MATLAB_COMMAND to the Matlab instance running NeuroElf
//
// http://SERVER:PORT/remote?gui=render_setview&guiargs=%27main%27%2C%5B15%2C45%2C0%2C0%2C1.2%5D
//
// passed on as
// neuroelf_gui('render_setview','main',[15,45,0,0,1.2]);

// settings and configurations

// path to Matlab binary
var MATLAB_PATH = '/Applications/MATLAB_R2013b.app/bin/matlab';
var MATLAB_ARGS = ['-nodesktop', '-nosplash', '-singleCompThread'];

// debugging
var NE_REMOTE_DEBUG = true;

// required modules
var http = require('http');
var spawn = require('child_process').spawn;
var url = require('url');

// spawn Matlab 
var matlab = spawn(MATLAB_PATH, MATLAB_ARGS);

// handle output of Matlab
matlab.stdout.on('data', function(matlab_output) {
    console.log(matlab_output.toString());
});

// log
console.log('MATLAB spawned.');

// test Matlab
matlab.stdin.write('a = 1; b = 2; c = a + b; disp(c);\n');

// start NeuroElf GUI in Matlab
matlab.stdin.write('neuroelf_gui;\n');

// create HTTP server
var http_server = http.createServer(function(req, res) {

    // debugging info
    if (NE_REMOTE_DEBUG) {
        console.log('HTTP request (' + req.method + '): ' + req.url);
    }

    // parse request
    var req_details = url.parse(req.url, true);

    // pass command on to Matlab?
    if (req_details.pathname === '/remote') {

        // get query
        var query = req_details.query;

        // has general command
        if (query.cmd) {

            // pass on command
            matlab.stdin.write(query.cmd + ';\n');

            // return status
            var response_body = 'command "' + query.cmd + '" passed on.\r\n';
            res.writeHead(200, {'Content-Length': response_body.length, 'Content-Type': 'text/plain'});
            res.write(response_body);
            res.end();
        }

        // has GUI command
        else if (query.gui) {

            // pass on GUI command
            var guicmd = 'neuroelf_gui(\'' + query.gui + '\'';
            if (query.guiargs)
                guicmd += ',' + query.guiargs;
            matlab.stdin.write(guicmd + ');\n');

            // return status
            var response_body = 'command "' + guicmd + ');" passed on.\r\n';
            res.writeHead(200, {'Content-Length': response_body.length, 'Content-Type': 'text/plain'});
            res.write(response_body);
            res.end();
        }

        // does not have command
        else {

            // return status
            var response_body = 'invalid request.';
            res.writeHead(501, {'Content-Length': response_body.length, 'Content-Type': 'text/plain'});
            res.write(response_body);
            res.end();
        }
    }

    // otherwise
    else {

        // return status
        var response_body = 'invalid request.';
        res.writeHead(501, {'Content-Length': response_body.length, 'Content-Type': 'text/plain'});
        res.write(response_body);
        res.end();
    }
});
http_server.listen(2280, '127.0.0.1');

