import re
import networkx as nx
import numpy as np
import os
import argparse

# MatrixMarket 파일을 읽고 값들을 교환한 후 정렬하는 함수
def swap_and_sort_values(data):
    # 값을 교환하고, 데이터를 오른쪽 값 기준으로 정렬
    swapped_data = [(u, v) for u, v in data]
    swapped_data.sort(key=lambda x: (x[0], x[1]))
    return swapped_data

# mtx 파일을 읽는 함수
def read_mtx_file(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
    header = [line for line in lines if line.startswith('%')]
    lines = [line for line in lines if not line.startswith('%')]
    
    size_info = lines.pop(0).split()
    num_nodes, _, num_edges = map(int, size_info)
    
    edges = []
    for line in lines:
        i, j = map(int, line.split())
        edges.append((i-1, j-1))
    
    return num_nodes, edges, header

# mtx 파일을 쓰는 함수
def write_mtx_file(graph, file_path, header):
    with open(file_path, 'w') as f:
        f.writelines(header)
        f.write(f'{graph.number_of_nodes()} {graph.number_of_nodes()} {graph.number_of_edges()}\n')
        
        # 그래프의 엣지를 리스트로 가져와서 값을 교환하고 정렬
        edges = list(graph.edges())
        sorted_edges = swap_and_sort_values(edges)
        #print(sorted_edges)
        
        for u, v in sorted_edges:
            f.write(f'{v + 1} {u + 1}\n')

def main(input_file):
    file_dir = f'./datasets/{input_file}.mtx'
    num_nodes, edges, header = read_mtx_file(file_dir)

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
        write_mtx_file(random_graph, output_file_path, header)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate random graphs using the configuration model.')
    parser.add_argument('input_file', type=str, help='Path to the input MTX file.')
    args = parser.parse_args()
    main(args.input_file)
