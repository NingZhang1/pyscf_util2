import re
from collections import namedtuple

# 定义数据结构
CQEDData = namedtuple("CQEDData", ["cmin", "ncsf", "ncfg", "roots"])
# 每个根的数据
RootData = namedtuple("RootData", ["eiCI", "ePT2", "etot"])


def parse_cqed_icipt2_output(text):
    """
    解析 cQED-iCIPT2 输出文本，提取所有块的信息。

    参数:
        text (str): 包含输出信息的字符串。

    返回:
        list of dict: 每个字典对应一个 Space 块，包含以下键：
            - spintwo (int): 自旋多重度
            - irrep (int): 不可约表示编号
            - data (list of CQEDData): 每个 CQEDData 包含 cmin, ncsf, ncfg 和 roots 列表，
                                         roots 列表中每个元素是 RootData (eiCI, ePT2, etot)
    """
    lines = text.splitlines()
    blocks = []
    current_block = None
    in_cqed = False

    # 正则表达式
    space_pattern = re.compile(r"spintwo\s+(\d+)\s+irrep\s+(\d+)")
    # 匹配数据行：以 | 开头，包含括号
    data_line_pattern = re.compile(
        r"^\|\s*([\d.eE+-]+)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*(.*)\s*\|$"
    )
    # 匹配括号内容：括号内三个浮点数，可能多个括号
    bracket_pattern = re.compile(
        r"\(\s*([-\d.eE+-]+)\s*,\s*([-\d.eE+-]+)\s*,\s*([-\d.eE+-]+)\s*\)"
    )

    for line in lines:
        line = line.strip()
        if not line:
            continue

        if "cQED-iCIPT2" in line:
            in_cqed = True
            continue

        if not in_cqed:
            continue

        # 匹配 Space 行，开始一个新块
        if "Space :" in line:
            match = space_pattern.search(line)
            if match:
                spintwo = int(match.group(1))
                irrep = int(match.group(2))
                # 如果已经有未关闭的块，先保存
                if current_block is not None:
                    blocks.append(current_block)
                current_block = {"spintwo": spintwo, "irrep": irrep, "data": []}
            continue

        # 跳过空行和分隔线
        if line.startswith("-"):
            continue

        # 匹配数据行
        match = data_line_pattern.match(line)
        if match:
            cmin = float(match.group(1))
            ncsf = int(match.group(2))
            ncfg = int(match.group(3))
            content = match.group(4)  # 括号部分，可能多个

            # 解析所有括号
            roots = []
            for bracket_match in bracket_pattern.finditer(content):
                eiCI = float(bracket_match.group(1))
                ePT2 = float(bracket_match.group(2))
                etot = float(bracket_match.group(3))
                roots.append(RootData(eiCI, ePT2, etot))

            if current_block is not None:
                current_block["data"].append(CQEDData(cmin, ncsf, ncfg, roots))
            continue

        # 遇到分隔线，可能结束当前块，但继续解析
        if line.startswith("---"):
            continue

    # 循环结束后保存最后一个块
    if current_block is not None:
        blocks.append(current_block)

    return blocks


def parse_cqed_icipt2_file(filepath, encoding="utf-8"):
    """
    从文件读取 cQED-iCIPT2 输出并解析。

    参数:
        filepath (str): 文件路径
        encoding (str): 文件编码，默认 utf-8

    返回:
        list of dict: 同 parse_cqed_icipt2_output，若失败返回 None
    """
    try:
        with open(filepath, "r", encoding=encoding) as f:
            content = f.read()
        return parse_cqed_icipt2_output(content)
    except FileNotFoundError:
        print(f"错误：文件 '{filepath}' 未找到。")
        return None
    except Exception as e:
        print(f"读取文件时发生错误：{e}")
        return None


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        filepath = sys.argv[1]
    else:
        filepath = "cqed_02S_0_0.out.dipole"

    result = parse_cqed_icipt2_file(filepath)
    if result is not None:
        # 简单打印结果示例
        for idx, block in enumerate(result):
            print(f"Block {idx}: spintwo={block['spintwo']}, irrep={block['irrep']}")
            for data in block["data"]:
                print(f"  cmin={data.cmin}, ncsf={data.ncsf}, ncfg={data.ncfg}")
                for i, root in enumerate(data.roots):
                    print(
                        f"    Root {i}: eiCI={root.eiCI:.12f}, ePT2={root.ePT2:.12f}, etot={root.etot:.12f}"
                    )
