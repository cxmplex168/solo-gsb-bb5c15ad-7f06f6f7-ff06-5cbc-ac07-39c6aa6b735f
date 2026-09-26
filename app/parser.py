"""化学式 / 物种解析器。

支持的语法（递归下降文法）::

    species    := formula [charge] | electron
    formula    := term+
    term       := '(' formula ')' [digits] | element [digits]
    element    := 大写字母 [小写字母]
    digits     := 非零开头的十进制正整数
    charge     := '^' [digits] ('+' | '-')
    electron   := 'e^-'

规则：

* 元素必须是真实存在的元素符号（118 种）；
* 下标必须为正整数，不允许前导 0；
* 圆括号最多嵌套四层；
* 重复出现的元素计数累加，括号乘以整体后同样累加；
* 电荷以后缀形式写在化学式末尾：``^+``、``^-``、``^n+``、``^n-``
  （n 为无前导零的正整数，省略时按 1 计）；无后缀表示中性；
* 电子仅写作 ``e^-``，其原子计数为空、电荷为 -1；
* 不支持点号结晶水（``CuSO4·5H2O``）与裸 ``+``/``-`` 电荷记号，
  方括号/花括号同样视为非法字符，必须改写为圆括号或拆成独立物质。
"""

from __future__ import annotations

MAX_PAREN_DEPTH = 4


class FormulaSyntaxError(ValueError):
    """非法化学式。

    ``position`` 为 0 起始的字符下标，指向首个非法字符
    （未闭合的 '(' 指向其本身，未配对的 ')' 指向 ')'，空串为 0；
    电荷后缀不完整时指向应有后续字符的位置，可等于字符串长度）。
    """

    def __init__(self, message: str, position: int):
        super().__init__(f"{message} (位置 {position})")
        self.message = message
        self.position = position


# 元素周期表 1–118 号元素符号
ELEMENT_SYMBOLS = frozenset(
    """
    H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni
    Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I
    Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt
    Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr
    Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og
    """.split()
)

# 只允许作为元素符号开头的大写字母范围（其余大写字母一定非法）
_UPPERCASE = frozenset("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
_LOWERCASE = frozenset("abcdefghijklmnopqrstuvwxyz")
_DIGITS = frozenset("0123456789")


def parse_species(text: str) -> tuple[dict[str, int], int]:
    """解析一个物种（化学式 + 可选电荷后缀，或电子 ``e^-``）。

    返回 ``({元素符号: 原子数}, 带符号电荷)``；电子的原子计数为空、
    电荷为 -1；无电荷后缀时电荷为 0。

    :raises FormulaSyntaxError: 化学式或电荷后缀非法，异常中带字符位置。
    """
    if not isinstance(text, str) or not text:
        raise FormulaSyntaxError("化学式不能为空", 0)

    # 电子：只允许整体写作 e^-（原子计数为空、电荷 -1）
    if text[0] == "e":
        if text == "e^-":
            return {}, -1
        if text.startswith("e^"):
            position = 2 if len(text) > 2 else len(text)
            if text.startswith("e^-") and len(text) > 3:
                position = 3  # e^- 之后的多余字符
            raise FormulaSyntaxError("电子仅写作 e^-", position)
        raise FormulaSyntaxError("元素符号必须以大写字母开头", 0)

    counts, pos = _parse_group(text, 0, depth=1, open_paren_pos=0)
    if pos < len(text) and text[pos] == "^":
        charge, pos = _parse_charge(text, pos)
        if pos != len(text):
            raise FormulaSyntaxError("电荷后缀后存在多余字符", pos)
        return dict(sorted(counts.items())), charge
    if pos != len(text):  # 理论上不会发生，留作防御
        raise FormulaSyntaxError("无法解析的残余字符", pos)
    return dict(sorted(counts.items())), 0


def parse_formula(text: str) -> dict[str, int]:
    """解析化学式，返回 ``{元素符号: 原子数}``。

    对于重复出现的元素（含括号展开）进行计数累加。
    末尾的电荷后缀（如 ``^2-``）会被接受但不计入返回值，
    需要电荷时请使用 :func:`parse_species`。

    :raises FormulaSyntaxError: 化学式非法，异常中带字符位置。
    """
    counts, _charge = parse_species(text)
    return counts


def _parse_group(
    text: str, pos: int, depth: int, open_paren_pos: int
) -> tuple[dict[str, int], int]:
    """解析一个括号组（或整个公式），返回 (计数表, 结束位置)。

    顶层的 ``^`` 会结束解析并把位置留给电荷后缀；括号内的 ``^`` 非法。
    """
    group: dict[str, int] = {}
    saw_term = False

    while pos < len(text):
        ch = text[pos]

        if ch == "(":
            # depth = 当前已开括号数 + 1；再开一层使已开数变为 depth，
            # depth > 4 即第 5 层，拒绝。
            if depth > MAX_PAREN_DEPTH:
                raise FormulaSyntaxError(
                    f"圆括号嵌套超过 {MAX_PAREN_DEPTH} 层", pos
                )
            inner, new_pos = _parse_group(text, pos + 1, depth + 1, pos)
            multiplier, new_pos = _parse_number(text, new_pos)
            _add(group, inner, multiplier)
            pos = new_pos
            saw_term = True

        elif ch == ")":
            if depth == 1:
                raise FormulaSyntaxError("存在未与 '(' 配对的 ')'", pos)
            if not saw_term:
                raise FormulaSyntaxError("括号内为空", pos)
            return group, pos + 1

        elif ch == "^":
            if depth > 1:
                raise FormulaSyntaxError(
                    "电荷后缀 '^' 只能出现在整个化学式的末尾", pos
                )
            if not saw_term:
                raise FormulaSyntaxError("电荷后缀 '^' 前缺少化学式", pos)
            break  # 顶层电荷后缀，交给 parse_species 处理

        elif ch in _UPPERCASE:
            symbol = ch
            pos += 1
            if pos < len(text) and text[pos] in _LOWERCASE:
                symbol += text[pos]
                pos += 1
            if symbol not in ELEMENT_SYMBOLS:
                raise FormulaSyntaxError(f"未知的元素符号 '{symbol}'", pos - len(symbol))
            number, pos = _parse_number(text, pos)
            group[symbol] = group.get(symbol, 0) + number
            saw_term = True

        elif ch in _LOWERCASE:
            raise FormulaSyntaxError("元素符号必须以大写字母开头", pos)

        elif ch in _DIGITS:
            raise FormulaSyntaxError("数字下标前缺少元素或右括号", pos)

        else:
            hint = _illegal_char_hint(ch)
            raise FormulaSyntaxError(hint, pos)

    if depth > 1:
        raise FormulaSyntaxError("存在未闭合的 '('", open_paren_pos)
    if not saw_term:
        raise FormulaSyntaxError("化学式不能为空", open_paren_pos)
    return group, pos


def _parse_number(text: str, pos: int) -> tuple[int, int]:
    """读取正整数下标；无下标时按 1 处理。"""
    if pos >= len(text) or text[pos] not in _DIGITS:
        return 1, pos
    if text[pos] == "0":
        raise FormulaSyntaxError("下标必须为正整数且不能以 0 开头", pos)
    start = pos
    while pos < len(text) and text[pos] in _DIGITS:
        pos += 1
    return int(text[start:pos]), pos


def _parse_charge(text: str, pos: int) -> tuple[int, int]:
    """解析 ``^`` 起始的电荷后缀，返回 (带符号电荷, 结束位置)。

    语法为 ``^`` 后接可选的电荷数（无前导零正整数，缺省为 1）与
    ``+``/``-`` 号；后缀之后是否还有字符由调用方检查。
    """
    pos += 1  # 跳过 '^'

    magnitude = 1
    if pos < len(text) and text[pos] in _DIGITS:
        if text[pos] == "0":
            raise FormulaSyntaxError("电荷数必须为正整数且不能以 0 开头", pos)
        start = pos
        while pos < len(text) and text[pos] in _DIGITS:
            pos += 1
        magnitude = int(text[start:pos])

    if pos >= len(text) or text[pos] not in ("+", "-"):
        raise FormulaSyntaxError("电荷后缀 '^' 后缺少 '+' 或 '-'", pos)
    sign = 1 if text[pos] == "+" else -1
    return sign * magnitude, pos + 1


def _add(target: dict[str, int], source: dict[str, int], multiplier: int) -> None:
    for symbol, amount in source.items():
        target[symbol] = target.get(symbol, 0) + amount * multiplier


def _illegal_char_hint(ch: str) -> str:
    if ch == ".":
        return "不支持点号结晶水，请将其拆成反应物/生成物中的独立物质"
    if ch in ("·", "•", "⋅"):
        return "不支持中点结晶水记号，请将其拆成反应物/生成物中的独立物质"
    if ch in ("+", "-"):
        return "裸 '+/-' 电荷记号不受支持，电荷请写作末尾的 ^+、^-、^n+、^n- 后缀"
    if ch in ("[", "]", "{", "}"):
        return f"不支持 '{ch}'，请改用圆括号 '(' ')'"
    if ch.isspace():
        return "化学式中不允许出现空白字符"
    if ord(ch) < 128:
        return f"化学式中存在非法字符 '{ch}'"
    return f"化学式中存在非 ASCII 字符 '{ch}'"
