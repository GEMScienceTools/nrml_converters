import sys
import os
import shutil
import runpy


def gem_rmtree(path, is_mandatory=False):
    try:
        shutil.rmtree(path)
    except OSError:
        if is_mandatory:
            raise
        else:
            pass


def gem_unlink(path, is_mandatory=False):
    try:
        os.unlink(path)
    except OSError:
        if is_mandatory:
            raise
        else:
            pass


def gem_run_script(prog, args):
    argv_orig = sys.argv
    sys.argv = [prog] + args

    syspath_orig = sys.path
    sys.path.insert(0, os.path.dirname(prog))

    try:
        runpy.run_path(prog, run_name='__main__')
    except SystemExit as excp:
        sys.argv = argv_orig
        sys.path = syspath_orig
        if excp.code == 0:
            pass
        else:
            raise
    sys.argv = argv_orig
    sys.path = syspath_orig
