import numpy as np
import galois

def bch_decode():
    # 參數
    m = 6
    n = 2**m - 1  # 63
    t = 2
    k = n - m*t   # 63 - 12 = 51

    simplified = 2  # 對應你程式裡 simplified != 1 的那條路

    # 建 GF(2^m)
    GF2 = galois.GF(2)
    prim_poly = galois.Poly.Int(67, field=GF2)   # 對應 MATLAB: D^6 + D + 1
    GF = galois.GF(2**m, irreducible_poly=prim_poly)

    # GF = galois.GF(2**m)

    # 常數
    alpha = GF(2)          # 對應 MATLAB: alpha = gf(2, m)
    zero  = GF(0)
    one   = GF(1)

    # input_string = "010011100001111110110100101010001110011010011010101110100000011"
    input_string = "010010100001111110110110101010001110011010111010101100100000011"
    if len(input_string) != n:
        raise ValueError(f"輸入字串長度 ({len(input_string)}) 與 BCH 碼長度 n={n} 不匹配。")

    # 字串 -> GF 向量
    reccode_v = np.array([int(ch) for ch in input_string], dtype=int)
    reccode   = GF(reccode_v)

    # 建 alpha_tb: [alpha^(2t), alpha^(2t-1), ..., alpha^1]
    alpha_tb = GF.Zeros(2*t)
    for i in range(2*t):
        alpha_tb[i] = alpha ** (2*t - i)

    # syndrome 計算 (跟你原本一樣用遞推)
    syndrome = GF.Zeros(2*t)
    for i in range(n):
        syndrome = syndrome * alpha_tb + reccode[i]

    # ===== Berlekamp–Massey =====
    # imba 初始化
    lambda_poly = GF.Zeros(t+1)
    lambda_poly[0] = one          # [1, 0, 0]
    b = GF.Zeros(t+2)
    b[1] = one                     # [0, 1, 0, 0]
    k_val = 0
    gamma = one
    delta = zero
    syndrome_array = GF.Zeros(t+1)

    if simplified == 1:
        # 你原本的 simplified branch，如果之後要用再翻
        for r in range(1, t+1):
            r1 = 2*t - 2*r + 2
            r2 = min(r1 + t, 2*t)
            # MATLAB 是 1-based, Python 0-based：要 -1
            r1_idx = r1 - 1
            r2_idx = r2          # slicing r1_idx:r2_idx 是 [r1_idx, r2_idx-1]
            num = r2 - r1 + 1

            syndrome_array[:num] = syndrome[r1_idx:r2_idx]
            delta = np.dot(syndrome_array, lambda_poly)

            lambda0 = lambda_poly.copy()
            # b2 那套我就先略，因為你現在 simplified=2
            raise NotImplementedError("simplified == 1 的路徑還沒翻完")
    else:
        # === 完整 BM 版本（你現在在用的 branch）===
        for r in range(1, 2*t + 1):
            # MATLAB:
            # r1 = 2*t - r + 1
            # r2 = min(r1 + t, 2*t)
            r1 = 2*t - r + 1
            r2 = min(r1 + t, 2*t)
            num = r2 - r1 + 1

            # 轉成 0-based index
            r1_idx = r1 - 1
            r2_idx = r2          # slicing 到 r2_idx-1

            selected_syndrome = syndrome[r1_idx:r2_idx]
            syndrome_array[:num] = selected_syndrome
            delta = np.dot(syndrome_array, lambda_poly)

            # 如果你還想 debug，可以 print 出來
            print(f"r = {r}, r1 = {r1}, r2 = {r2}, selected_syndrome = {syndrome_array}")
            print(f"lambda (before) = {lambda_poly}")
            print(f"b (before)      = {b}")
            print(f"delta           = {int(delta)}")
            print(f"gamma           = {int(gamma)}")
            print(f"k               = {k_val}")

            lambda0 = lambda_poly.copy()
            lambda_poly = gamma * lambda_poly - delta * b[0:t+1]

            if (delta != zero) and (k_val >= 0):
                b[1:1+t+1] = lambda0  # b(2:2+t) = lambda0 (t+1 elements)
                gamma = delta
                k_val = -k_val - 1
            else:
                b[1:1+t+1] = b[0:t+1]
                # gamma 不變
                k_val = k_val + 1

            print(f"lambda (after)  = {lambda_poly}")
            print(f"b (after)       = {b}")
            print(f"gamma (after)   = {int(gamma)}")
            print(f"k (after)       = {k_val}")
            print("="*40)

    # ===== Chien search =====

    # inverse_tb(i) = alpha^(-i+1), i=1..t+1
    inverse_tb = GF.Zeros(t+1)
    for i in range(1, t+1+1):
        inverse_tb[i-1] = alpha ** (-i + 1)

    # accu_tb 初始為全 1
    accu_tb = GF.Ones(t+1)
    error = np.zeros(n, dtype=int)

    for i in range(n):
        lambda_v = np.dot(lambda_poly, accu_tb)
        if lambda_v == zero:
            error[i] = 1
        else:
            error[i] = 0
        # element-wise 乘 inverse_tb
        accu_tb = accu_tb * inverse_tb

    found = np.nonzero(error != 0)[0]  # index from 0..n-1
    print("Error positions (0-based) =", found)

    return found, error

if __name__ == "__main__":
    found, error = bch_decode()

