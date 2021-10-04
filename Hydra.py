import numpy as np
import matplotlib.pyplot as plt

version = 2
error_correction_level = 'L'
module = 25
mask = '000'
V = np.zeros((module, module), dtype=int)


class GF():  # 域内的乘除运算与初始化
    alpha_to_int = {}  # 对数表示转十进制
    int_to_alpha = {}  # 十进制转对数表示

    def __init__(self):
        num = 1
        for idx in range(255):
            self.alpha_to_int[idx] = num
            self.int_to_alpha[num] = idx
            num = num * 2
            if num > 255:
                num = num ^ 285

    def mul(self, a, b):
        if a == 0 or b == 0:
            return 0
        return self.alpha_to_int[(self.int_to_alpha[a] + self.int_to_alpha[b])
                                 % 255]

    def div(self, a, b):
        if a == 0:
            return 0
        x = (self.int_to_alpha[a] - self.int_to_alpha[b]) % 255
        if x < 0:
            x += 255
        return self.alpha_to_int[x]


gf = GF()


class polynomial():  #time：多项式的次数 Poly：将多项式由低次向高次储存 mul/mod：计算该多项式与另一多项式的乘/除法
    time = 0
    Poly = []

    def to_alpha(self):  # 全部转化为对数表示
        for i in range(self.time + 1):
            self.Poly[i] = gf.int_to_alpha[self.Poly[i]]

    def to_int(self):  # 全部转化为十进制表示
        for i in range(self.time + 1):
            self.Poly[i] = gf.alpha_to_int[self.Poly[i]]

    def mul(self, oth):
        ans = polynomial()
        ans.Poly = []
        ans.time = self.time + oth.time
        for i in range(ans.time + 1):
            temp = 0
            to = min(self.time, i)
            for j in range(to + 1):
                if i - j > oth.time:
                    continue
                temp ^= gf.mul(self.Poly[j], oth.Poly[i - j])
            ans.Poly.append(temp)
        return ans

    def mod(self, oth):  # 计算除法，需先转化为倒序储存，返回余数
        if self.time < oth.time:
            return self
        for i in range(self.time - oth.time + 1):
            temp = gf.div(self.Poly[i], oth.Poly[0])
            if temp == 0:
                continue
            for j in range(oth.time + 1):
                oth.Poly[j] = gf.mul(oth.Poly[j], temp)
                self.Poly[i + j] ^= oth.Poly[j]
        self.Poly = self.Poly[self.time - oth.time + 1:]
        self.time = oth.time - 1
        return self


def input_version_and_error_correction_level(Version, Error_correction_level):
    global version
    version = Version
    global error_correction_level
    error_correction_level = Error_correction_level
    global module
    module = 17 + 4 * version
    global V
    V = np.zeros((module, module), dtype=int)


def char_to_val(x):  # 使用Alphanumeric编码
    dict = {
        ' ': 36,
        '$': 37,
        '%': 38,
        '*': 39,
        '+': 40,
        '-': 41,
        '.': 42,
        '/': 43,
        ':': 44
    }
    if '0' <= x <= '9':
        return int(x)
    elif 'A' <= x <= 'Z':
        return ord(x) - ord('A') + 10
    else:
        return dict[x]


def get_number_of_data_bits():
    dict = {
        '1L': 152,
        '1M': 128,
        '1Q': 104,
        '1H': 72,
        '2L': 272,
        '2M': 224,
        '2Q': 176,
        '2H': 128
    }
    return dict[str(version) + error_correction_level]


def get_number_of_ECC():
    dict = {
        '1L': 7,
        '1M': 10,
        '1Q': 13,
        '1H': 17,
        '2L': 10,
        '2M': 16,
        '2Q': 22,
        '2H': 28
    }
    return dict[str(version) + error_correction_level]


def coding_bit_stream(str):  # 使用Alphanumeric获得的数据流的头部信息和数据部分
    ans = '0010'  # 编码方式
    length = len(str)
    ans = ans + '{:09b}'.format(length)  # 长度转二进制，不足9位补0
    for idx in range(0, length, 2):
        if idx == length - 1:
            continue
        ans = ans + '{:011b}'.format(  # 数据串中两位字符一组转二进制，不足11位补0
            char_to_val(str[idx]) * 45 + char_to_val(str[idx + 1]))
    if length & 1:  # 多出来一个字符
        ans = ans + '{:06b}'.format(char_to_val(str[length - 1]))
    return ans


def get_bit_stream(str):  #获得数据流
    words = coding_bit_stream(str)
    data_bits = get_number_of_data_bits()
    length = len(words)
    if data_bits - length >= 4:  # 终止符
        words = words + '0000'
        length += 4
    words = words + (8 - length % 8) * '0'  # 八位1组补足
    length += (8 - length % 8)
    while length + 16 <= data_bits:  # 补足码字数
        words = words + '1110110000010001'
        length += 16
    if (length < data_bits):
        words = words + '11101100'
    return words


def get_message_polynomial(bit_stream):  #获得信息多项式
    Mx = polynomial()
    Mx.Poly = []
    length = len(bit_stream)
    Mx.time = (length // 8 - 1)
    for idx in range(0, length, 8):
        Mx.Poly.append(int(bit_stream[idx:idx + 8], 2))  # 八位一个码字
    Mx.Poly = Mx.Poly[::-1]  # 低次到高次排列
    return Mx


def get_generator_polynomial():  #获得生成多项式
    m = get_number_of_ECC()
    Gx = polynomial()
    Gx.time = 1
    Gx.Poly = [1, 1]
    for i in range(1, m):
        temp = polynomial()
        temp.time = 1
        temp.Poly = [gf.alpha_to_int[i], 1]
        Gx = Gx.mul(temp)
    return Gx


def get_error_correcting_code(bit_stream):  #获得纠错码
    Mx = get_message_polynomial(bit_stream)
    Gx = get_generator_polynomial()
    m = get_number_of_ECC()
    Mx.Poly = Mx.Poly[::-1]
    Gx.Poly = Gx.Poly[::-1]
    for i in range(m):
        Mx.Poly.append(0)  # 对齐
    Mx.time += m
    Mx = Mx.mod(Gx)
    words = ''
    for i in Mx.Poly:
        words = words + '{:08b}'.format(i)
    return words


def combine_BS_ECC(bit_stream,
                   error_correcting_code):  # 将数据流码字和纠错码码字转换为最终的01数据串
    words = bit_stream + error_correcting_code
    if (version == 2):
        words += '0' * 7
    return words


def draw_signal():  #填充 QR 码定位标志、校正标志、定时标志等非数据位置

    #绘制定位标志
    for i in range(7):
        V[i][0] = V[i][6] = V[i][module - 1] = V[i][module - 7] = 1
        V[module - 1 - i][0] = V[module - 1 - i][6] = 1
        V[0][i] = V[6][i] = V[module - 1][i] = V[module - 7][i] = 1
        V[0][module - 1 - i] = V[6][module - 1 - i] = 1
    for i in range(3):
        for j in range(3):
            V[2 + i][2 + j] = 1
            V[module - 5 + i][2 + j] = 1
            V[2 + i][module - 5 + j] = 1

    #绘制校正标志
    if version == 2:
        for i in range(5):
            V[module - 9][module - 9 + i] = 1
            V[module - 9 + i][module - 9] = 1
            V[module - 5 - i][module - 5] = 1
            V[module - 5][module - 5 - i] = 1
        V[module - 7][module - 7] = 1
    #绘制定时标志
    for i in range(8, module - 8, 2):
        V[6][i] = V[i][6] = 1


def draw_data(data):  #填充数据串
    length = len(data)
    cnt = 0
    p = -1
    x = module
    y = module - 1
    while cnt < length:
        if (version == 2 and x == module - 4 and y == module - 9):  #校正标志左侧特殊一列
            for i in range(5):
                x -= 1
                V[x][y - 1] = data[cnt]
                cnt += 1
            continue
        if ((x == module - 1 and p == 1) or (x == 0 and p == -1)  #触碰边界或定位标志
                or (x == 9 and (y < 8 or y >= module - 8) and p == -1)
                or (x == module - 9 and y <= 8 and p == 1)):
            y -= 2
            p *= -1
            if x == module - 1 and y == 8:
                x -= 8
        elif (version == 2 and (x == module - 4 and y == module - 5)):  #触碰校正标志
            x -= 6
        elif (version == 2 and (x == module - 10 and y == module - 7)):
            x += 6
        elif (x + p == 6):  #触碰计时标志
            x += p * 2
        elif (x == 9 and y == 8):
            y -= 3
            p *= -1
        else:  #正常情况
            x += p
        V[x][y] = data[cnt]
        if cnt == length - 1:
            break
        V[x][y - 1] = data[cnt + 1]
        cnt += 2


def count_mask(x, y, mask):
    if mask == 0:
        return int((x + y) % 2 == 0)
    elif mask == 1:
        return int((x % 2 == 0))
    elif mask == 2:
        return int(y % 3 == 0)
    elif mask == 3:
        return int((x + y) % 3 == 0)
    elif mask == 4:
        return int((x // 2 + y // 3) % 2 == 0)
    elif mask == 5:
        return int(((x * y) % 2 + (x * y) % 3) == 0)
    elif mask == 6:
        return int(((x * y) % 2 + (x * y) % 3) % 2 == 0)
    else:
        return int(((x + y) % 2 + (x * y) % 3) % 2 == 0)


def apply_mask(mask):  #应用掩码
    if type(mask) == str:
        mask = int(mask, 2)
    p = -1
    x = module
    y = module - 1
    while True:
        if (version == 2 and x == module - 4 and y == module - 9):  #校正标志左侧特殊一列
            for i in range(5):
                x -= 1
                V[x][y - 1] ^= count_mask(x, y - 1, mask)
            continue
        if ((x == module - 1 and p == 1) or (x == 0 and p == -1)  #触碰边界或定位标志
                or (x == 9 and (y < 8 or y >= module - 8) and p == -1)
                or (x == module - 9 and y <= 8 and p == 1)):
            y -= 2
            p *= -1
            if x == module - 1 and y == 8:
                x -= 8
        elif (version == 2 and (x == module - 4 and y == module - 5)):  #触碰校正标志
            x -= 6
        elif (version == 2 and (x == module - 10 and y == module - 7)):
            x += 6
        elif (x + p == 6):  #触碰计时标志
            x += p * 2
        elif (x == 9 and y == 8):
            y -= 3
            p *= -1
        else:
            x += p
        V[x][y] ^= count_mask(x, y, mask)
        V[x][y - 1] ^= count_mask(x, y - 1, mask)
        if x == module - 9 and y == 1:
            break


def evaluate_mask():  #评估掩码
    tot = 0
    #rule 1 连续同色色条检测
    for i in range(module):
        tmp = 0
        lst = -1
        for j in range(module):
            if V[i][j] == lst:
                tmp += 1
            else:
                lst = V[i][j]
                if tmp >= 5:
                    tot += tmp - 4
                tmp = 0
        if tmp >= 5:
            tot += tmp - 4
            tmp = 0

    for i in range(module):
        tmp = 0
        lst = -1
        for j in range(module):
            if V[j][i] == lst:
                tmp += 1
            else:
                lst = V[j][i]
                if tmp >= 5:
                    tot += tmp - 4
                tmp = 0
        if tmp >= 5:
            tot += tmp - 4
            tmp = 0

    #rule 2 同色色块检测
    for i in range(module - 1):
        for j in range(module - 1):
            if (V[i][j] == V[i][j + 1] == V[i + 1][j] == V[i + 1][j + 1]):
                tot += 3

    #rule 3 特殊模式图案检测 10111010000
    for i in range(module):
        for j in range(module):
            tmp = ''
            if (j < module - 10):
                for k in range(11):
                    tmp += str(V[i][j + k])
                if (int(tmp, 2) == 1488):
                    tot += 40
            tmp = ''
            if (j >= 10):
                for k in range(11):
                    tmp += str(V[i][j - k])
                if (int(tmp, 2) == 1488):
                    tot += 40
            tmp = ''
            if (i < module - 10):
                for k in range(11):
                    tmp += str(V[i + k][j])
                if (int(tmp, 2) == 1488):
                    tot += 40
            tmp = ''
            if (i >= 10):
                for k in range(11):
                    tmp += str(V[i - k][j])
                if (int(tmp, 2) == 1488):
                    tot += 40

    #rule 4 黑白比例检测
    cnt = 0
    for i in range(module):
        for j in range(module):
            cnt += V[i][j]
    percntage = cnt * 100 // (module**2)
    tot += 10 * min(
        abs(50 - percntage // 5 * 5) * 2,
        abs(50 - (percntage // 5 + 1) * 5) * 2)
    return tot


def get_mask():  #选择最佳掩码
    mask_id = 0
    mask_evaluation = -1
    for i in range(8):
        apply_mask(i)
        evaluation = evaluate_mask()
        apply_mask(i)
        if (evaluation < mask_evaluation or mask_evaluation == -1):
            mask_evaluation = evaluation
            mask_id = i
    return '{:03b}'.format(mask_id)


def get_format_inf(mask):  #获得格式信息串
    level_code = {'L': '01', 'M': '00', 'Q': '11', 'H': '10'}
    format_inf = level_code[error_correction_level] + mask
    Mx = []
    for i in format_inf:
        Mx.append(int(i))
    for i in range(10):
        Mx.append(0)
    Gx = [1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1]
    length = len(Mx)
    for i in range(length - 10):  # 求余
        if Mx[i] != 0:
            for j in range(11):
                Mx[i + j] ^= Gx[j]
    for i in Mx[length - 10:]:
        format_inf = format_inf + str(i)
    return '{:015b}'.format(int(format_inf, 2) ^ int('101010000010010', 2))


def draw_format_inf(format_inf):  #填充格式信息
    for i in range(6):
        V[8][i] = V[module - 1 - i][8] = format_inf[i]
    V[8][7] = V[module - 7][8] = format_inf[6]
    V[8][8] = V[8][module - 8] = format_inf[7]
    V[7][8] = V[8][module - 7] = format_inf[8]
    for i in range(6):
        V[5 - i][8] = V[8][module - 6 + i] = format_inf[9 + i]
    V[module - 8][8] = 1


def get_QR_code(Str, version, error_correction_level):  #获得一个二维码的二维数组
    input_version_and_error_correction_level(version, error_correction_level)
    bit = get_bit_stream(Str)
    ecc = get_error_correcting_code(bit)
    ans = combine_BS_ECC(bit, ecc)
    draw_signal()
    draw_data(ans)
    mask = get_mask()
    apply_mask(mask)
    format_inf = get_format_inf(mask)
    draw_format_inf(format_inf)
    return V


def draw_QR_code(Str, version, error_correction_level):  #画一个二维码
    plt.axis('off')
    plt.imshow(get_QR_code(Str, version, error_correction_level),
               cmap='gray_r')
    plt.show()
