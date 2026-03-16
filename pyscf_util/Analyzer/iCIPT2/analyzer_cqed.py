from pyscf_util.Analyzer.iCIPT2.analyzer import *


def parse_cqed_icipt2_output(text):
    """
    解析 cQED-iCIPT2 输出文本，提取 spintwo, irrep 和表格数据。

    参数:
        text (str): 包含输出信息的字符串。

    返回:
        dict: 包含以下键的字典：
            - spintwo (int): 自旋多重度
            - irrep (int): 不可约表示编号
            - data (list of tuple): 每个元组包含
              (cmin, ncsf, ncfg, eiCI, ePT2, etot)
              其中 cmin 为浮点数，ncsf/ncfg 为整数，其余为浮点数。
    """

    lines = text.splitlines()

    spintwo = None
    irrep = None
    data = []

    # 正则表达式匹配 "Space : spintwo 0 irrep 0" 中的数字
    space_pattern = re.compile(r"spintwo\s+(\d+)\s+irrep\s+(\d+)")

    find_cqed = False

    for line in lines:
        line = line.strip()

        if "cQED-iCIPT2" in line:
            find_cqed = True

        if not find_cqed:
            continue

        # 提取 spintwo 和 irrep
        if "Space :" in line:
            match = space_pattern.search(line)
            if match:
                spintwo = int(match.group(1))
                irrep = int(match.group(2))
            continue

        # 跳过空行和分隔线（以 --- 开头）
        if not line or line.startswith("-"):
            continue

        # 数据行以 '|' 开头且包含 '('
        if line.startswith("|") and "(" in line:
            parts = line.split("|")
            # parts[0] 为空，parts[1]~parts[3] 为前三列，parts[4] 为括号部分，parts[5] 为空
            if len(parts) >= 5:
                cmin_str = parts[1].strip()
                ncsf_str = parts[2].strip()
                ncfg_str = parts[3].strip()
                paren_str = parts[4].strip()

                # 转换前三个数值

                try:
                    cmin = float(cmin_str)
                    ncsf = int(ncsf_str)
                    ncfg = int(ncfg_str)
                except Exception as e:
                    continue

                # 处理括号内的三个数：移除首尾括号，按逗号分割
                inner = paren_str.strip("()").split(",")
                if len(inner) == 3:
                    eiCI = float(inner[0].strip())
                    ePT2 = float(inner[1].strip())
                    etot = float(inner[2].strip())

                    # data.append((cmin, ncsf, ncfg, eiCI, ePT2, etot))

                    data.append(
                        iCIPT2_Data(
                            ncfg=ncfg,
                            ncsf=ncsf,
                            evar=eiCI,
                            ept=ePT2,
                            etot=etot,
                        )
                    )

    return {"spintwo": spintwo, "irrep": irrep, "data": data}


def parse_cqed_icipt2_file(filepath, encoding="utf-8"):
    """
    从文件读取 cQED-iCIPT2 输出并解析。

    参数:
        filepath (str): 文件路径
        encoding (str): 文件编码，默认 utf-8

    返回:
        dict: 解析结果，同 parse_cipt2_output
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

    print(parse_cqed_icipt2_file("cqed_cc-pvtz_3_1_quadrupole.out"))
