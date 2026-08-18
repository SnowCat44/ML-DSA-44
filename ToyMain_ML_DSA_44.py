import importlib.util
import os
import sys

MODULE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Toy_ML_DSA_44.py")


def load_dsa_module(path=MODULE_PATH):
    spec = importlib.util.spec_from_file_location("toy_ml_dsa", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    dsa = load_dsa_module()
    msg = b"hello world"

    # 1) KeyGen의 pk, sk 정상 반환 여부 확인
    pk, sk = dsa.ML_DSA_KeyGen()
    print(f"[1] KeyGen -> pk: {'OK' if pk is not None else 'FAIL'}, sk: {'OK' if sk is not None else 'FAIL'}")

    # 2) Sign의 sigma 정상 반환 여부 확인
    sigma = dsa.ML_DSA_Sign(sk, msg)
    print(f"[2] Sign   -> sigma: {'OK' if sigma is not None else 'FAIL'}")

    # 3) Verify의 정상 동작 여부 확인
    ok = dsa.ML_DSA_Verify(pk, msg, sigma)
    print(f"[3] Verify -> {'PASS' if ok else 'FAIL'}")

    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
