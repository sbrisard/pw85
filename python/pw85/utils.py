import configparser
import ctypes
import os.path
import sys

c_double_p = ctypes.POINTER(ctypes.c_double)


def get_config_dir():
    home = os.path.expanduser("~")
    if sys.platform.startswith("darwin"):
        return os.path.join(home, "Library", "Application Support", "pw85")
    elif sys.platform.startswith("win32"):
        return os.path.join(os.environ["APPDATA"], "pw85")
    elif sys.platform.startswith("linux"):
        return os.path.join(home, ".pw85")


def get_config_option(name):
    path = os.path.join(get_config_dir(), "pw85.cfg")
    if os.path.isfile(path):
        cfg = configparser.ConfigParser()
        cfg.read(path)
        return cfg["pw85"][name]
    else:
        raise RuntimeError("Cannot file configuration file: {}".format(path))


def load_library(name):
    return ctypes.cdll.LoadLibrary(get_config_option(name))
