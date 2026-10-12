"""运行所有单元测试.

无需 pytest, 用 assert 即可. 退出码: 0=全过, 1=有失败.
"""
import os
import sys
import traceback

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _run(mod):
    print(f'\n=== {mod} ===')
    try:
        __import__(mod)
        print(f'  PASS')
        return True
    except SystemExit as e:
        if e.code == 0:
            print(f'  PASS')
            return True
        print(f'  FAIL: SystemExit {e.code}')
        return False
    except Exception as e:
        print(f'  FAIL: {type(e).__name__}: {e}')
        traceback.print_exc()
        return False


def main():
    failed = 0
    for mod in ['tests.test_io_bruker', 'tests.test_io_sdf', 'tests.test_retrieve', 'tests.test_score', 'tests.test_bridge_dept', 'tests.test_new_features', 'tests.test_io_bruker_2d', 'tests.test_correlation_2d', 'tests.test_algorithmic_assignment']:
        if not _run(mod):
            failed += 1
    print()
    if failed:
        print(f'FAILED: {failed} module(s)')
        sys.exit(1)
    else:
        print('ALL PASSED')
        sys.exit(0)


if __name__ == '__main__':
    main()