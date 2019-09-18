'''Checks the accuracy of various implementations of the function f.'''
import os.path
import shutil

import h5py
import numpy as np
import requests

import pypw85

PROXIES = {'http': 'http://proxy.enpc.fr:3128',
           'https': 'https://proxy.enpc.fr:3128'}

REF_DATA_URL = 'https://zenodo.org/record/3323683/files/pw85_ref_data-20190712.h5'
REF_DATA_PATH = 'pw85_ref_data.h5'


def download_ref_data():
    if not os.path.exists(REF_DATA_PATH):
        r = requests.get(REF_DATA_URL, proxies=PROXIES, stream=True)
        if r.status_code == 200:
            with open(REF_DATA_PATH, 'wb') as f:
                shutil.copyfileobj(r.raw, f)
        else:
            raise RuntimeError('could not retrieve reference data')


def accuracy_histogram(func, spheroids, directions, lambdas, expecteds):
    prec = np.finfo(np.float64).precision
    hist = np.zeros(prec+2, dtype=np.uint64)
    edges = np.linspace(-prec-1., 1., num=prec+3)
    for i1, q1 in enumerate(spheroids):
        print("{}/{}".format(i1+1, spheroids.shape[0]))
        for i2, q2 in enumerate(spheroids):
            for i, r12_i in enumerate(directions):
                for j, lambda_j in enumerate(lambdas):
                        exp = expecteds[i1, i2, i, j]
                        if exp == 0.:
                            raise RuntimeError('expected value should not be zero')
                        act = func(lambda_j, r12_i, q1, q2)
                        if np.isnan(act):
                            print('nan')
                            act = 1e100
                        err = np.abs((act-exp)/exp)

                        if err == 0.:
                            hist[0] += 1
                        else:
                            bin_index = int(np.floor(np.log10(err)+prec+1))
                            if bin_index < 0:
                                bin_index = 0
                            if bin_index > prec+1:
                                bin_index = prec+1
                            hist[bin_index] += 1
    return hist, edges


if __name__ == '__main__':
    download_ref_data()
    with h5py.File(REF_DATA_PATH, 'r') as f:
        spheroids = np.array(f['spheroids'])
        directions = np.array(f['directions'])
        lambdas = np.array(f['lambdas'])
        expecteds = np.array(f['F'])

    for func, path in [(pypw85.f, 'implementation-01.csv'),
                       (pypw85.f2, 'implementation-02.csv')]:

        hist, edges = accuracy_histogram(func, spheroids, directions, lambdas, expecteds)
        centers = 0.5*(edges[1:]+edges[:-1])
        widths = edges[1:]-edges[:-1]
        with open(path, 'w') as f:
            for x, dx, y in zip(centers, widths, hist):
                f.write('{},{},{}\n'.format(x, dx, y))
