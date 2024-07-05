import subprocess
import re
import os
import numpy as np

def create_edges_file(mtx_file_path, edges_file_path):
    """mtx 파일을 .edges 파일로 변환하여 노드 숫자를 -1씩 감소 및 정렬"""
    with open(mtx_file_path, 'r') as file:
        lines = file.readlines()

    # 첫 번째 줄을 "노드 수 엣지 수" 형식으로 변경
    node_edge_info = lines[1].strip().split()
    if len(node_edge_info) < 3:
        print(f"Error: Unexpected format in line 2 of {mtx_file_path}")
        return
    
    node_count = node_edge_info[0]
    edge_count = node_edge_info[2]

    with open(edges_file_path, 'w') as file:
        file.write(f"{node_count} {edge_count}\n")
        # 노드 숫자를 -1씩 감소시키고 작은 숫자가 왼쪽에 오도록 정렬
        for i in range(2, len(lines)):
            nodes = list(map(int, lines[i].strip().split()))
            if len(nodes) == 2:
                nodes = [node - 1 for node in nodes]
                nodes.sort()
                file.write(f"{nodes[0]} {nodes[1]}\n")

def run_subgraph_counts(file_path):
    """subgraph_counts.py를 실행하고 결과를 반환"""
    command = ['python', './escape/wrappers/subgraph_counts.py', file_path, '5', '-i']
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1)

    output = ""
    for line in process.stdout:
        output += line
        # 진행 상황 출력
        if "progress" in line.lower():
            print(line.strip())
    
    process.wait()
    if process.returncode != 0:
        print(f"Error running subgraph_counts.py: {process.stderr.read()}")
    return output


def parse_output(output):
    """subgraph_counts.py의 출력을 파싱하여 매핑된 패턴과 값을 딕셔너리로 반환"""
    mapping = {
        "4-1": "4-clique",
        "4-2": "Diamond",
        "4-3": "Tailed triangle",
        "4-4": "4-cycle",
        "4-5": "3-star",
        "4-6": "3-path",
        "5-1": "4-star",
        "5-2": "Prong",
        "5-3": "4-path",
        "5-4": "Forktailed-tri",
        "5-5": "Lontailed-tri",
        "5-6": "Doubletailed-tri",
        "5-7": "Tailed-4-cycle",
        "5-8": "5-cycle",
        "5-9": "Hourglass",
        "5-10": "Cobra",
        "5-11": "Stingray",
        "5-12": "Hatted-4-cycle",
        "5-13": "3-wedge-col",
        "5-14": "3-tri-collision",
        "5-15": "Tailed-4-clique",
        "5-16": "Triangle-strip",
        "5-17": "Diamond-wed-col",
        "5-18": "4-wheel",
        "5-19": "Hatted-4-clique",
        "5-20": "Almost-5-clique",
        "5-21": "5-clique"
    }
    
    reverse_mapping = {v: k for k, v in mapping.items()}
    pattern_counts = {}

    for line in output.split("\n"):
        match = re.match(r"(\w[\w\s\-]*)\s+(\d+)\s+(\d+)\s+([\d\.]+)", line)
        if match:
            pattern_name = match.group(1).strip()
            count = float(match.group(3).strip())
            if pattern_name in reverse_mapping:
                mapped_name = reverse_mapping[pattern_name]
                if mapped_name not in pattern_counts:
                    pattern_counts[mapped_name] = []
                pattern_counts[mapped_name].append(count)

    return pattern_counts

def average_counts(pattern_counts):
    """각 패턴의 값을 평균내어 반환"""
    averaged_counts = {pattern: np.mean(counts) for pattern, counts in pattern_counts.items()}
    return averaged_counts

def process_files(base_file_path, random_folder_path):
    # .edges 파일 생성
    base_edges_file_path = f'{base_file_path}.edges'
    create_edges_file(base_file_path, base_edges_file_path)

    # Base 파일 처리
    print("Base File processing...")
    output = run_subgraph_counts(base_edges_file_path)
    base_counts = parse_output(output)

    # Base 파일 결과 출력
    print("Base File Results:")
    print("{")
    for key in sorted(base_counts.keys()):
        print(f'    "{key}" : {base_counts[key][0]:.2f},')  # base file has single value for each pattern
    print("}")

    # .edges 파일 삭제
    if os.path.exists(base_edges_file_path):
        os.remove(base_edges_file_path)

    # Random graph files 처리
    random_counts = {}
    for i in range(1, 11):
        random_mtx_file_path = os.path.join(random_folder_path, f'{os.path.basename(base_file_path).split(".")[0]}_r{i}.mtx')
        random_edges_file_path = f'{random_mtx_file_path}.edges'
        create_edges_file(random_mtx_file_path, random_edges_file_path)
        output = run_subgraph_counts(random_edges_file_path)
        counts = parse_output(output)

        # .edges 파일 삭제
        if os.path.exists(random_edges_file_path):
            os.remove(random_edges_file_path)

        # 각 패턴의 값을 합산
        for pattern, count_list in counts.items():
            if pattern not in random_counts:
                random_counts[pattern] = []
            random_counts[pattern].extend(count_list)

    # 평균 계산
    averaged_counts = average_counts(random_counts)

    # Random graph 파일 결과 출력
    print("Random Graphs Results (Averaged):")
    print("{")
    for key in sorted(averaged_counts.keys()):
        print(f'    "{key}" : {averaged_counts[key]:.2f},')
    print("}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python main.py <base_filename>")
        sys.exit(1)

    base_filename = sys.argv[1]
    base_file_path = f'./datasets/{base_filename}.mtx'
    random_folder_path = f'./datasets/{base_filename}_random'

    process_files(base_file_path, random_folder_path)
