import importlib.util
import os
import sys

MODULE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "민규_Toy_ML_DSA.py")


def load_dsa_module(path=MODULE_PATH):
    spec = importlib.util.spec_from_file_location("toy_ml_dsa", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    dsa = load_dsa_module()

    pk, sk = dsa.ML_DSA_KeyGen()
    msg = b"hello world"
    sig = dsa.ML_DSA_Sign(sk, msg)
    ok = dsa.ML_DSA_Verify(pk, msg, sig)

    print(f"메시지: {msg!r}")
    print(f"서명 검증 결과: {'PASS' if ok else 'FAIL'}")

    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
