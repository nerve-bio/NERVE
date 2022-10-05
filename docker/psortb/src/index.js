const fs = require('fs');
const {spawn} = require("child_process");
var Mutex = require('async-mutex').Mutex;

const express = require("express");
const app = express();
app.use(express.json({limit: '50mb'}));

const mutex = new Mutex();
let id = 0;

function encodeParameters(parameters) {
    let re = '';
    for (let parameter in parameters) {
        if (parameters[parameter] !== false) {
            re += ' -' + parameter;
            if (parameters[parameter] !== '' && parameters[parameter] !== true) {
                re += ' ' + parameters[parameter];
            }
        }
    }
    return re;
}

app.post('/execute', async (req, res) => {
    if (
        req.body.gram && ['p', 'n', 'a'].indexOf(req.body.gram)>=0
        && req.body.seq
    ) {
        let uuid = 0;
        await mutex.runExclusive(async () => {
            uuid = id++;
        });
        const workingDirectory = __dirname + '/temp_wd/' + uuid;
    	const inputFilePath = workingDirectory + '/input.fasta'
    	const resultFilePath = workingDirectory + '/result.txt'
        
    	const gram = req.body.gram;
    	const args = {
    	    p: gram==='p',
    	    n: gram==='n',
    	    a: gram==='a',
    	    c: !isNaN(req.body.cutoff)?req.body.cutoff:false,
    	    d: !isNaN(req.body.divergent)?req.body.divergent:false,
    	    f: req.body.format?req.body.format:false,
    	    e: req.body.exact?req.body.exact:false,
    	    o: req.body.output?req.body.output:false,
    	    v: true,
    	    l: workingDirectory,
    	    "-seq": inputFilePath
    	};
    
    	
    	// init working directory
	fs.mkdir(workingDirectory, {recursive: true}, err => {
	    if (err) {
	        console.error(err);
	        res.status(500).send('error while creating working directory');
	    } else {
    	        // put seq file
    	        fs.writeFile(inputFilePath, req.body.seq, err => {
	            if (err) {
	                console.error(err);
	                res.status(500).send('error while creating input file');
	            } else {
	                // start child psortb process
	                const command = '/usr/local/psortb/bin/psort' + encodeParameters(args);
                        console.log('Spawning child process with command: ' + command);
    
                        const commandSplit = command.split(' ');
	                const child = spawn(commandSplit[0], commandSplit.slice(1), {cwd: workingDirectory});
	                let stdout = '';
	                let stderr = '';
                        child.stdout.on('data', (data) => {
                            stdout += data;
                            console.log("child stdout:" + data);
                        });
                        child.stderr.on('data', (data) => {
                            stderr += data;
                            console.error("child stderr:" + data);
                        });
                        child.on('exit', async (code, signal) => {
        	            // read results
        	            fs.readFile(resultFilePath, 'utf8', (err, result) => {
                                if (err) {
                                    console.error(err);
	                            res.status(500).send('error while reading result file');
                                } else {
                                    // cleanup
                                    fs.rm(workingDirectory, {recursive: true, force: true}, err => {
                                        if (err) {
                                            console.error(err);
	                                    res.status(500).send('error while cleaning working directory');
                                        } else {
                                            // send results
                                            res.send({
                                                stdout: stdout,
                                                stderr: stderr,
                                                result: result
                                            });
                                        }
                                    });
                                }
                             });
                        });
                        child.on('error', (err) => {
                            res.status(500).send('error while spawning psortb child process');
                            console.error(err);
                        });
	            }
	        });
	    }
	});
    } else {
        res.status(400).send('bad request');
    }
});

app.listen(8080);
console.log('Listening on port 8080');
