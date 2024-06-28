import networkx as nx
import numpy as np
import os
import argparse

def read_mtx_file(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
    header = [line for line in lines if line.startswith('%')]
    lines = [line for line in lines if not line.startswith('%')]
    
    size_info = lines.pop(0).split()
    num_rows, num_cols, num_entries = map(int, size_info)
    
    edges = []
    for line in lines:
        i, j = map(int, line.split())
        edges.append((i-1, j-1))
    
    return num_rows, edges

def write_mtx_file(graph, file_path):
    with open(file_path, 'w') as f:
        f.write('%%MatrixMarket matrix coordinate pattern symmetric\n')
        f.write(f'{graph.number_of_nodes()} {graph.number_of_nodes()} {graph.number_of_edges()}\n')
        
        for u, v in graph.edges():
            f.write(f'{u + 1} {v + 1}\n')

def main(input_file):
    file_dir = f'./datasets/{input_file}.mtx'
    num_nodes, edges = read_mtx_file(file_dir)

    G = nx.Graph()
    G.add_edges_from(edges)

    # 차수 분포 추출
    degree_sequence = [d for n, d in G.degree()]

    # 랜덤 그래프를 저장할 폴더 생성
    output_dir = f'./datasets/{input_file}_random'
    os.makedirs(output_dir, exist_ok=True)

    # 랜덤 그래프 10개 생성 및 저장
    for i in range(10):
        idx = i + 1
        random_graph = nx.configuration_model(degree_sequence, create_using=nx.Graph())
        random_graph = nx.Graph(random_graph)
        random_graph.remove_edges_from(nx.selfloop_edges(random_graph))

        output_file_path = os.path.join(output_dir, f'{input_file}_r{idx}.mtx')
        write_mtx_file(random_graph, output_file_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate random graphs using the configuration model.')
    parser.add_argument('input_file', type=str, help='Path to the input MTX file.')
    args = parser.parse_args()
    main(args.input_file)
