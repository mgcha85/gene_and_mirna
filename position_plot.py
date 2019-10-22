import socket
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np


class position_plot:
    def __init__(self):
        hostname = socket.gethostname()
        print(hostname)
        if hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinfomatics'
        elif hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        elif 'evc' in hostname:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        else:
            print('wrong option')
            return
        self.cell_lines = ['K562', 'HepG2', 'A549', 'GM12878', 'HEK293']

    def region_plot(self, xstart, xend, ystart, yend, color='b', label='tag'):
        width = xend - xstart
        height = yend - ystart
        return matplotlib.patches.Rectangle((xstart, ystart), width, height, alpha=0.8, color=color, label=label)
    
    def point_plot(self, ax, point, color='r', label='tss'):
        ax.plot([point, point], [0, 1.5], linewidth=2, alpha=0.5, color=color, label=label)
        
    def run(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5', 'Search_TSSs.xlsx')
        df = pd.read_excel(fpath)

        for idx in df.index:
            print('{:,d} / {:,d}'.format(idx + 1, df.shape[0]))
            fig = plt.figure()
            ax = fig.add_subplot(111)

            cell_line = df.loc[idx, 'cell_line']
            gene = df.loc[idx, 'gene']
            strand = df.loc[idx, 'strand']
            tss = df.loc[idx, 'tss']
            tag_starts = map(int, df.loc[idx, 'tag_starts'].split(';'))
            tag_ends = map(int, df.loc[idx, 'tag_ends'].split(';'))
            gene_tss = df.loc[idx, 'gene_tss']

            ttl = 'gene: [{}],  cell_line: {}, strand: {}'.format(gene, cell_line, strand)
            starts = []
            for i, (lstart, lend, label, color) in enumerate(zip(['pre_start', 'gene_start'], ['pre_end', 'gene_end'], ['pre-miRNA', 'gene'], ['g', 'c'])):
                start = df.loc[idx, lstart]
                end = df.loc[idx, lend]
                starts.append(start)

                rect = self.region_plot(start, end, 0, 1, color=color, label=label)
                ax.add_patch(rect)
                for j, (tstart, tend) in enumerate(zip(tag_starts, tag_ends)):
                    starts.append(tstart)
                    if j == 0:
                        label = 'tag'
                    else:
                        label = None
                    rect = self.region_plot(tstart, tend, 0, 0.5, label=label)
                    ax.add_patch(rect)

            self.point_plot(ax, tss, color='r', label='pre_miRNA_tss')
            self.point_plot(ax, gene_tss, color='k', label='gene_tss')
            if strand == '+':
                range_start = min(starts) - 500
                range_end = df.loc[idx, 'pre_end'] + 500
            else:
                range_start = df.loc[idx, 'pre_start'] - 500
                range_end = max(starts) + 500

            length = abs(gene_tss - tss)
            label_xpos = np.mean([gene_tss, tss]).astype(int)
            ax.text(label_xpos, 1.2, "{:,d}bp".format(length), ha='center', size=8)
            ax.plot(sorted([gene_tss, tss]), [1.4, 1.4], color='k')
            plt.xlim([range_start, range_end])
            plt.ylim([0, 5])
            plt.legend()
            plt.title(ttl)
            plt.savefig('figures/{}_{}.png'.format(gene, cell_line))
            plt.close()


if __name__ == '__main__':
    pp = position_plot()
    pp.run()
