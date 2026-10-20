"""
DCM_fusion_main_preopgrouped.py - runs DCM_fusion_main.py unchanged, except that patients are grouped by
their PreOp category (CATEGORY_SOURCE = 'PreOp', applied in memory - DCM_fusion_main.py itself is not edited).
DCM_fusion_main.py writes to fusion_vs_lami_results/category_preop/; every file it produces there is then
renamed with '_preopgrouped' appended to its name.

Run: python DCM_fusion_main_preopgrouped.py
"""

import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
MAIN = os.path.join(HERE, 'DCM_fusion_main.py')
SUFFIX = '_preopgrouped'

src = open(MAIN, encoding='utf-8').read()
src, n = re.subn(r"^CATEGORY_SOURCE = 'PostOp'", "CATEGORY_SOURCE = 'PreOp'", src, flags=re.M)
assert n == 1, "expected exactly one CATEGORY_SOURCE assignment, found {}".format(n)
exec(compile(src, MAIN, 'exec'), {'__name__': '__main__', '__file__': MAIN})

out_dir = os.path.join(HERE, 'fusion_vs_lami_results', 'category_preop')
renamed = 0
for root, _, files in os.walk(out_dir):
    for f in files:
        stem, ext = os.path.splitext(f)
        if stem.endswith(SUFFIX):      # already suffixed (e.g. the explicit preop-grouped difference plot)
            continue
        os.replace(os.path.join(root, f), os.path.join(root, stem + SUFFIX + ext))
        renamed += 1
print("\nAppended '{}' to {} file(s) in {}".format(SUFFIX, renamed, out_dir))
