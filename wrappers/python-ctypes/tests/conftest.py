import pytest


def pytest_addoption(parser):
    parser.addini("http_proxy", "")
    parser.addini("https_proxy", "")
