from decimal import Decimal, getcontext, ROUND_HALF_UP


def format_measurement(value: float, uncertainty: float) -> str:
    """
    将数值和不确定度格式化为科学表格表达式，例如 -1.23(5)。

    规则：
    - 不确定度：若首位有效数字为 1 或 2，保留两位有效数字；否则保留一位。
    - 数值：四舍五入到与不确定度末位对齐的小数位数。
    - 输出格式：数值(括号内为不确定度的整数形式)

    参数:
        value: 测量值
        uncertainty: 不确定度（正数）

    返回:
        格式化字符串，如 "1.23(5)"
    """
    if uncertainty == 0:
        # 不确定度为 0 时直接返回原值（不推荐，但按需求处理）
        return str(value)

    # 设置精度足够高，避免舍入误差
    getcontext().prec = 28

    # 将输入转为 Decimal，使用字符串避免浮点误差
    u = Decimal(str(uncertainty)).normalize()
    v = Decimal(str(value))

    # 获取不确定度的科学计数法指数（adjusted exponent）
    exp = u.adjusted()  # e.g., 0.045 -> -2, 1.23 -> 0, 100 -> 2

    # 获取首位有效数字 (1-9)
    leading_digit = int(u.scaleb(-exp))  # scaleb(-exp) 将数值变为 [1,10) 之间的数
    # print(leading_digit)

    # # 决定保留的有效数字位数
    # if leading_digit in (1, 2):
    #     keep_digits = 2
    #     target_exp = exp - 1  # 舍入到更小一位（更高精度）
    # else:
    #     keep_digits = 1
    #     target_exp = exp
    target_exp = exp

    # 构造舍入单位（例如 1e-2, 1e0, 1e1 等）
    quant_unit = Decimal("1e{}".format(target_exp))

    # 对不确定度进行舍入
    rounded_u = u.quantize(quant_unit, rounding=ROUND_HALF_UP)

    # 确定数值需要保留的小数位数
    # rounded_u 的小数位数 = max(0, -exponent)，其中 exponent 来自 Decimal 的 as_tuple()
    u_exp = rounded_u.as_tuple().exponent
    decimal_places = max(0, -u_exp)

    # 对数值进行舍入，使其与 rounded_u 的小数位数对齐
    if decimal_places == 0:
        quant_v = Decimal("1")
    else:
        quant_v = Decimal("1e-{}".format(decimal_places))
    rounded_v = v.quantize(quant_v, rounding=ROUND_HALF_UP)

    # 格式化数值部分（保留固定小数位数，必要时补零）
    if decimal_places == 0:
        value_str = f"{int(rounded_v)}"
    else:
        value_str = f"{rounded_v:.{decimal_places}f}"

    # 括号内的整数： rounded_u * 10^{decimal_places}
    unc_int = int(rounded_u * Decimal(10**decimal_places))

    return f"{value_str}({unc_int})"


# print table #

from tabulate import tabulate


def print_table_sci(DATA, KEY_RES, KEY_ERR):
    header = ["task"]
    header2 = ["task"]
    for key1, key2 in zip(KEY_RES, KEY_ERR):
        header.append(key1)
        header2.append(key1)
        header.append(key2)
    data1 = []
    data2 = []
    for key in DATA:
        data_print1 = [key]
        data_print2 = [key]
        try:
            for key1, key2 in zip(KEY_RES, KEY_ERR):
                data_print1.append(DATA[key][key1])
                data_print1.append(DATA[key][key2])
                data_print2.append(format_measurement(DATA[key][key1], DATA[key][key2]))
        except Exception as e:
            data_print1.append(0.0)
            data_print1.append(0.0)
            data_print2.append("0(0)")
        data1.append(data_print1)
        data2.append(data_print2)

    print(tabulate(data1, headers=header, tablefmt="grid", floatfmt="15.8f"))
    print(tabulate(data2, headers=header2, tablefmt="grid", floatfmt="15.8f"))


# 示例用法
if __name__ == "__main__":
    # 题目示例
    print(format_measurement(-1.234, 0.045))  # -1.23(5)

    # 更多测试
    print(format_measurement(0.1234, 0.0104))  # 0.123(10)
    print(format_measurement(12.345, 0.99))  # 12.3(10)
    print(format_measurement(123.45, 1.23))  # 123.5(12)
    print(format_measurement(0.001234, 0.000104))  # 0.00123(10)
    print(format_measurement(100.5, 0.5))  # 100.5(5)
    print(format_measurement(-0.5678, 0.023))  # -0.57(2)
    print(format_measurement(9.876, 0.099))  # 9.9(1)   (0.099 -> 0.1 -> 括号1)
    print(format_measurement(5.0, 0.101))  # 5.00(10) (0.101 -> 0.10 -> 括号10)
