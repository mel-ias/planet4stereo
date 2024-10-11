import os, sys, subprocess

def run_cmd(bin, args, **kw):
    """
    wrapper around subprocess function to excute bash commands
    Parameters
    ----------
    bin: str
        command to be excuted (e.g., stereo or gdalwarp)
    args: list
        arguments to the command as a list
    Retuns
    ----------
    out: str
        log (stdout) as str if the command executed, error message if the command failed
    """
    #Note, need to add full executable
    #from dshean/vmap.py
    #binpath = os.path.join('/home/sbhushan/src/StereoPipeline/bin',bin)
    binpath = find_executable(bin)
    if binpath is None:
        msg = ("Unable to find executable %s\n" 
        #"Install ASP and ensure it is in your PATH env variable\n" "https://ti.arc.nasa.gov/tech/asr/intelligent-robotics/ngt/stereo/" 
        % bin)
        sys.exit(msg)
    #binpath = os.path.join('/opt/StereoPipeline/bin/',bin)
    call = [binpath,]
    call.extend(args)
    #print("my call", call)
    #print(' '.join(call))
    try:
        out = subprocess.run(call,encoding='UTF-8') #Mel: changed subprocess.check_call to the newer subprocess.run (this waits for the result and pipes all printing stuff to terminal)
    except:
        out = "the command {} failed to run, see corresponding asp log".format(call)
    return out

def find_executable(executable, path=None):
    """Tries to find 'executable' in the directories listed in 'path'.

    A string listing directories separated by 'os.pathsep'; defaults to
    os.environ['PATH'].  Returns the complete filename or None if not found.
    """
    _, ext = os.path.splitext(executable)
    if (sys.platform == 'win32') and (ext != '.exe'):
        executable = executable + '.exe'

    if os.path.isfile(executable):
        return executable

    if path is None:
        path = os.environ.get('PATH', None)
        if path is None:
            try:
                path = os.confstr("CS_PATH")
            except (AttributeError, ValueError):
                # os.confstr() or CS_PATH is not available
                path = os.defpath
    if not path:
        return None

    paths = path.split(os.pathsep)
    for p in paths:
        f = os.path.join(p, executable)
        if os.path.isfile(f):
            # the file exists, we have a shot at spawn working
            return f
    return None
