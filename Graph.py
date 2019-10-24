import networkx as nx
import pandas as pd
import pickle as pkl
import socket


class Graph:
    def __init__(self):
        hostname = socket.gethostname()
        if hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinformatics'
        elif hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'

    def run(self):
        df = pd.read_excel('correlation_report2.xlsx')
        df_grp = df.groupby('GENE')
        G = nx.Graph()

        for gene, df_sub in df_grp:
            print(gene)
            G.add_node(gene)

            for idx in df_sub.index:
                mirna = df_sub.loc[idx, 'miRNA']
                G.add_node(mirna)
                weigth = df_sub.loc[idx, 'value']
                G.add_edge(gene, mirna, weight=weigth)

        with open('graph.cha', 'wb') as f:
            pkl.dump(G, f)

    def plot(self):
        import matplotlib.pyplot as plt
        with open('graph.cha', 'rb') as f:
            G = pkl.load(f)

        nx.draw(G, with_labels=True)
        nx.draw_shell(G, with_labels=True)
        plt.show()


if __name__ == '__main__':
    g = Graph()
    g.run()
