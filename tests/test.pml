load ref.center.pdb
load cmp.center.pdb

pair_fit cmp.center & not H*, ref.center & not H*

python
from pymol import *
print cmd.get_object_matrix("cmp.center")

python end
