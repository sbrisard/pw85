import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    name = 'histograms.csv'
    data = np.loadtxt(name, delimiter=',')
    bin_width = 1.0
    bar_width = 0.4*bin_width

    fig, ax = plt.subplots(figsize=(6., 4.))
    ax.bar(data[:, 0]-0.5*bar_width, data[:, 1], bar_width, label='Implementation 1')
    ax.bar(data[:, 0]+0.5*bar_width, data[:, 2], bar_width, label='Implementation 2')
    ax.legend()

    ax.set_xticks(data[:, 0])
    xticklabels = list(data[:, 0].astype(np.int))
    xticklabels[0] = '≤ 0'
    xticklabels[-1] = '≥ {}'.format(data.shape[0]-1)
    ax.set_xticklabels(xticklabels, rotation=45)

    ax.set_xlabel('Precision (number of digits)')
    ax.set_ylabel('Count')

    fig.tight_layout()
    fig.savefig('histograms.pdf')
    fig.savefig('histograms.png', dpi=600)
