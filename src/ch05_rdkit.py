# %%

from rdkit.Chem import Draw
from rdkit.Chem import EnumerateHeterocycles
from rdkit.Chem import AllChem
import itertools
import copy
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, TemplateAlign, rdDepictor, rdFMCS
from rdkit.Chem.Draw import IPythonConsole

IPythonConsole.ipython_useSVG = True
rdDepictor.SetPreferCoordGen(True)

# %% 構造を描画してみよう

mol = Chem.MolFromSmiles("c1ccccc1")
mol

# %%A 複数の化合物を一度に取り扱うには？

mols = Chem.SDMolSupplier("./data/ch05_compounds.sdf")
len(mols)

Draw.MolsToGridImage([mol for mol in mols])

# %% ヘテロシャッフリングをしてみる - 画像表示

sildenafil = Chem.MolFromSmiles(
    'CCCC1=NN(C)C2=C1NC(=NC2=O)C1=C(OCC)C=CC(=C1)S(=O)(=O)N1CCN(C)CC1')
vardenafil = Chem.MolFromSmiles(
    'CCCC1=NC(C)=C2N1NC(=NC2=O)C1=C(OCC)C=CC(=C1)S(=O)(=O)N1CCN(CC)CC1')
rdDepictor.Compute2DCoords(sildenafil)
rdDepictor.Compute2DCoords(vardenafil)
res = rdFMCS.FindMCS([sildenafil, vardenafil], completeRingsOnly=True,
                     atomCompare=rdFMCS.AtomCompare.CompareAny)
MCS = Chem.MolFromSmarts(res.smartsString)
rdDepictor.Compute2DCoords(MCS)

TemplateAlign.AlignMolToTemplate2D(sildenafil, MCS)
TemplateAlign.AlignMolToTemplate2D(vardenafil, MCS)
Draw.MolsToGridImage([sildenafil, vardenafil], legends=[
                     'sildenafil', 'vardenafil'])

# %% ヘテロシャッフリングをしてみる - 使ってみる


def check_mol(mol):
    arom_atoms = mol.GetAromaticAtoms()
    symbols = [atom.GetSymbol()
               for atom in arom_atoms if not atom.IsInRingSize(5)]
    if not symbols:
        return True
    elif 'O' in symbols or 'S' in symbols:
        return False
    else:
        return True


class HeteroShuffle():

    def __init__(self, mol, query):
        self.mol = mol
        self.query = query
        self.subs = Chem.ReplaceCore(self.mol, self.query)
        self.core = Chem.ReplaceSidechains(self.mol, self.query)
        self.target_atomic_nums = [6, 7, 8, 16]

    def make_connectors(self):
        n = len(Chem.MolToSmiles(self.subs).split('.'))
        map_no = n+1
        self.rxn_dict = {}
        for i in range(n):
            self.rxn_dict[i+1] = AllChem.ReactionFromSmarts(
                '[{0}*][*:{1}].[{0}*][*:{2}]>>[*:{1}][*:{2}]'.format(i+1, map_no, map_no+1))
        return self.rxn_dict

    def re_construct_mol(self, core):
        '''
        re construct mols from given substructures and core
        '''
        keys = self.rxn_dict.keys()
        ps = [[core]]
        for key in keys:
            ps = self.rxn_dict[key].RunReactants([ps[0][0], self.subs])
        mol = ps[0][0]
        try:
            smi = Chem.MolToSmiles(mol)
            mol = Chem.MolFromSmiles(smi)
            Chem.SanitizeMol(mol)
            return mol
        except:
            return None

    def get_target_atoms(self):
        '''
        get target atoms for replace
        target atoms means atoms which don't have anyatom(*) in neighbors
        '''
        atoms = []
        for atom in self.core.GetAromaticAtoms():
            neighbors = [a.GetSymbol() for a in atom.GetNeighbors()]
            if '*' not in neighbors and atom.GetSymbol() != '*':
                atoms.append(atom)
        print(len(atoms))
        return atoms

    def generate_mols(self):
        atoms = self.get_target_atoms()
        idxs = [atom.GetIdx() for atom in atoms]
        combinations = itertools.product(
            self.target_atomic_nums, repeat=len(idxs))
        smiles_set = set()
        self.make_connectors()
        for combination in combinations:
            target = copy.deepcopy(self.core)
            for i, idx in enumerate(idxs):
                target.GetAtomWithIdx(idx).SetAtomicNum(combination[i])
            smi = Chem.MolToSmiles(target)
            target = Chem.MolFromSmiles(smi)
            if target is not None:
                n_attachment = len(
                    [atom for atom in target.GetAtoms() if atom.GetAtomicNum() == 0])
                n_aromatic_atoms = len(list(target.GetAromaticAtoms()))
                if target.GetNumAtoms() - n_attachment == n_aromatic_atoms:
                    try:
                        mol = self.re_construct_mol(target)
                        if check_mol(mol):
                            smiles_set.add(Chem.MolToSmiles(mol))
                    except:
                        pass
        mols = [Chem.MolFromSmiles(smi) for smi in smiles_set]
        return mols


# Gefitinib
mol1 = Chem.MolFromSmiles(
    'COC1=C(C=C2C(=C1)N=CN=C2NC3=CC(=C(C=C3)F)Cl)OCCCN4CCOCC4')
core1 = Chem.MolFromSmiles('c1ccc2c(c1)cncn2')
# Oxaprozin
mol2 = Chem.MolFromSmiles('OC(=O)CCC1=NC(=C(O1)C1=CC=CC=C1)C1=CC=CC=C1')
core2 = Chem.MolFromSmiles('c1cnco1')

Draw.MolsToGridImage([mol1, mol2])

# %%

ht = HeteroShuffle(mol1, core1)
res = ht.generate_mols()
print(len(res))
Draw.MolsToGridImage(res, molsPerRow=5)

# %%

ht = HeteroShuffle(mol2, core2)
res = ht.generate_mols()
print(len(res))
Draw.MolsToGridImage(res, molsPerRow=5)

# %% EnumerateHeterocycles

omeprazole = Chem.MolFromSmiles(
    'CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C(N2)C=C(C=C3)OC')
enumerated_mols = EnumerateHeterocycles.EnumerateHeterocycles(omeprazole)
enumerated_mols = [m for m in enumerated_mols]
print(len(enumerated_mols))
# 384
Draw.MolsToGridImage(enumerated_mols[:6], molsPerRow=3)

# %%

ringinfo = omeprazole.GetRingInfo()
ringinfo.AtomRings()
# ((1, 6, 5, 4, 3, 2), (14, 13, 17, 16, 15), (18, 19, 20, 21, 15, 16))
protected_omeprazole = copy.deepcopy(omeprazole)
for idx in ringinfo.AtomRings()[1]:
    atom = protected_omeprazole.GetAtomWithIdx(idx)
    atom.SetProp("_protected", "1")
for idx in ringinfo.AtomRings()[2]:
    atom = protected_omeprazole.GetAtomWithIdx(idx)
    atom.SetProp("_protected", "1")
enumerated_mols2 = EnumerateHeterocycles.EnumerateHeterocycles(
    protected_omeprazole)
enumerated_mols2 = [m for m in enumerated_mols2]
print(len(enumerated_mols2))
# 4
Draw.MolsToGridImage(enumerated_mols2)

# %%
