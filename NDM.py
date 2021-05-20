import numpy as np
from sklearn.decomposition import NMF
import pandas as pd
import networkx as nx


def comput_nmf(v,n):
    model = NMF(n_components=n, init=None, random_state=0)
    W = model.fit_transform(v)
    H = model.components_
    return W,H


def initScore():
    S_File = r"subcellular_score.txt"
    H_File = r"homologous_score.txt"
    S_Score = dict()
    S_Value = []
    H_Score = dict()
    H_Value = []
    N_Score = dict()

    dataS = [line.strip('\n').split() for line in open(S_File, 'r')]
    for key, value in dataS:
        value = float(value)
        S_Value.append(value)
        S_Score[key] = value

    dataH = [line.strip('\n').split() for line in open(H_File, 'r')]
    for key, value in dataH:
        value = float(value)
        H_Value.append(value)
        H_Score[key] = value
    max_H_Value = max(H_Value)
    for key in H_Score.keys():
        H_Score[key] = H_Score[key] / max_H_Value
        # N_Score[key] = (H_Score[key] + S_Score[key]) / 2
        N_Score[key] = H_Score[key] * S_Score[key]
    return N_Score


def LoadEssential():
    filename1 = r"essential.txt"
    essentialProtein = []
    with open(filename1, 'r') as f:
        for line in f:
            # result.append(list(line.strip('\n').split(',')))
            # essetialProtein.append(line.strip('\n'))
            essentialProtein += line.strip('\n').split(' ')
    return essentialProtein


def AccEssential(score, essentialProtein, size):
    count = 0
    for index, value in score[:size]:
        if index in essentialProtein:
            count += 1
    
    return count

def saveFile(new_Score):
    
    score_array = [new_Score]
    save = pd.DataFrame(score_array)
    save.to_csv('test2.csv', index=False, header=True)  

def loadEdge(filename):
    G = nx.Graph()
    edges = [line.strip('\n').split() for line in open(filename, 'r')]
    # G = nx.read_weighted_edgelist(filename, create_using=nx.Graph)
    G.add_weighted_edges_from(edges)
    return G

def comScore(n, nodes):
    
    n = 850  
    file = ['co_neighbor.txt', 'Complex_weight.txt', 'Domain_weight.txt']
    tempData = listData = np.zeros((len(nodes), len(nodes)))
    for filename in file:
        subG = loadEdge(filename)
        for node in nodes:
            if node not in list(subG.nodes):
                subG.add_node(node)
        tempData = nx.to_numpy_matrix(subG, nodelist=nodes)
        listData += (tempData * 1 / 3)  
    v = listData
    W, H = comput_nmf(v, n) 
    new_v = np.dot(W, H)  


    new_edges = []   
    N = len(nodes) 
    sum_v = new_v.sum(axis=0)
    for i in range(N):
        for j in range(i, N):
            if sum_v[i] != 0:
                new_v[i][j] = max(new_v[i][j], new_v[j][i]) / sum_v[i]
                new_v[j][i] = new_v[i][j]
            else:
                new_v[j][i] = new_v[i][j] = 0  

    
    new_G = nx.Graph()
    new_G = nx.from_numpy_matrix(new_v)

    N_Score = dict()
    N_Score = initScore()
    for i in range(N):
        N_Score[i] = N_Score.pop(nodes[i])

    # alph = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    alph = [0.3]
    score_array = []

    for alp in alph:
        print("特征数为：%d,alpha为:%.1f" % (n, alp))
        new_Score = dict()
        new_Score = nx.pagerank(new_G, alpha=alp, nstart=N_Score, personalization=N_Score)
        for i in range(N):
            new_Score[nodes[i]] = new_Score.pop(i)
        # saveFile(new_Score)
        score_array.append(new_Score)
        sort_Score = sorted(new_Score.items(), key=lambda x: x[1], reverse=True)
        essentialProtein = LoadEssential()
        size = [100, 200, 300, 400, 500, 600]
        for si in size:
            acc_Count = AccEssential(sort_Score, essentialProtein, si)
            print("When size is %d, ACC num is %d" % (si, acc_Count))
    save = pd.DataFrame(score_array)
    save.to_csv('0209850result.csv', index=True, header=True)  # index=False,heade


def readCsv():
    G = nx.Graph()
    edges = [line.strip('\n').split() for line in open(r"DIP2010.txt", 'r')]
    G.add_edges_from(edges)
    nodes = list(G.nodes)
    filename1 = r"test2.csv"
    df = pd.read_csv(filename1, sep=',', header=None)
    dict_tmp = dict(zip(df.values[0, :], df.values[1, :]))
    # dict_tmp = dict(zip(nodes, df.values[1, :]))
    for item in dict_tmp.items():
        print(item)
    sort_Score = sorted(dict_tmp.items(), key=lambda x: x[1], reverse=True)
    essentialProtein = LoadEssential()
    size = [93, 200, 300, 400, 500, 600]
    for si in size:
        acc_Count = AccEssential(sort_Score, essentialProtein, si, nodes)
        print("When size is %d, ACC num is %d" % (si, acc_Count))


if __name__ == "__main__":
    G = nx.Graph()
    edges = [line.strip('\n').split() for line in open(r"DIP2010.txt", 'r')]
    G.add_edges_from(edges)
    nodes = list(G.nodes)

    n_featue = [850]
    score_array = []
    for n in n_featue:
        comScore(n, nodes)
