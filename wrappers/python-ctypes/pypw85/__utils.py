import configparser
import ctypes
import os.path
import sys

c_double_p = ctypes.POINTER(ctypes.c_double)


def pw85_data_dir():
    home = os.path.expanduser('~')
    if sys.platform.startswith('darwin'):
        return os.path.join(home, 'Library', 'Application Support', 'pw85')
    elif sys.platform.startswith('win32'):
        return os.path.join(os.environ['APPDATA'], 'pw85')
    elif sys.platform.startswith('linux'):
        return os.path.join(home, '.pw85')


def load_library(name):
    path = os.path.join(pw85_data_dir(), 'pw85.cfg')
    if os.path.isfile(path):
        cfg = configparser.ConfigParser()
        cfg.read(path)
        return ctypes.cdll.LoadLibrary(cfg["pw85"][name])
    else:
        raise RuntimeError("Cannot file configuration file: {}".format(path))
