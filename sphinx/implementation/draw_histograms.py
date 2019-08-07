import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    name = 'implementation-{:02d}.csv'
    data1, data2 = [np.loadtxt(name.format(i), delimiter=',') for i in range(1, 3)]
    # Assume constant bin width
    bin_width = data1[0, 1]
    bar_width = 0.4*bin_width

    fig, ax = plt.subplots(figsize=(6., 4.))
    for i, data, in enumerate((data1, data2), 1):
        ax.bar(data[:, 0]+(-1)**i*0.5*bar_width, data[:, 2], bar_width, label='Implementation #{}'.format(i))
    ax.legend()

    ax.set_xticks(data1[:, 0])
    xticklabels = list((data1[:, 0]+0.5*data1[:, 1]).astype(np.int))
    xticklabels.insert(0, '< {}'.format(xticklabels[0]))
    ax.set_xticklabels(xticklabels, rotation=45)

    ax.set_xlabel('Precision (number of digits)')
    ax.set_ylabel('Count')

    fig.tight_layout()
    fig.savefig('histograms.pdf')
    fig.savefig('histograms.png', dpi=600)
