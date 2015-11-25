import runpy
import sys


def gem_run_script(prog, args):
    argv_orig = sys.argv
    sys.argv = [prog] + args
    try:
        runpy.run_path(prog, run_name='__main__')
    except SystemExit as excp:
        sys.argv = argv_orig
        if excp.code == 0:
            pass
        else:
            raise
    sys.argv = argv_orig
