from __future__ import print_function
import sys, os, shutil, subprocess

def run(command, options=None):
    if isinstance(command, list):
        commandstring = ' '.join( map(lambda x: str(x), command))
        commandlist = command
    else:
        commandstring = command
        commandlist = command.split()
    # print("command=%s\ncommandstring=%s\ncommandlist=%s" %(command, commandstring, commandlist))
    if options==None:
        return subprocess.check_call(commandstring, shell=True)
    elif options=="file":
        try:
            stdout_file = open('./simfem.stdout', 'a+')
            stderr_file = open('./simfem.stderr', 'a+')
        except:
            raise IOError('cannot open files')
        proc = subprocess.Popen(args=commandlist, stdout=stdout_file, stderr=stderr_file, shell=False)
        returncode = proc.wait()
        if returncode:
            raise RuntimeError("error in command '%s'" %commandstring)
        return returncode
    else:
        p = subprocess.Popen(commandlist, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout, stderr = p.communicate()
        if p.returncode != 0:
            raise RuntimeError("error in command '%s' (running in %s) errorcode=%d\n\nstdout=\x1b[6;30;42m %s \x1b[0m" %(commandstring, os.getcwd(), p.returncode,stdout.decode('ascii')))
        return p.returncode
