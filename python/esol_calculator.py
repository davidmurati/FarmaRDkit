"""
ESOL Calculator - Estimación de Solubilidad Acuosa
Implementación de la ecuación de Delaney (ESOL):
  logS = 0.16 - 0.63*cLogP - 0.0062*MW + 0.066*RB - 0.74*AP
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import math


SOLUBILITY_CLASSES = [
    (-12, -10, "Insoluble",       "#ef4444", "< 0.01 μg/mL"),
    (-10, -6,  "Poco soluble",    "#f97316", "0.01–10 μg/mL"),
    (-6,  -4,  "Moderadamente soluble", "#eab308", "10–100 μg/mL"),
    (-4,  -2,  "Soluble",         "#22c55e", "0.1–10 mg/mL"),
    (-2,   0,  "Muy soluble",     "#3b82f6", "10–1000 mg/mL"),
    (0,    14, "Altamente soluble","#8b5cf6", "> 1000 mg/mL"),
]


def classify_solubility(logs: float) -> dict:
    """Clasifica la solubilidad basada en logS."""
    for low, high, label, color, desc in SOLUBILITY_CLASSES:
        if low <= logs < high:
            return {"class": label, "color": color, "description": desc}
    return {"class": "Indefinida", "color": "#6b7280", "description": "Fuera de rango"}


def logs_to_mg_per_l(logs: float, mol_weight: float) -> float:
    """Convierte logS (mol/L) a mg/L."""
    try:
        mol_per_l = 10 ** logs
        return mol_per_l * mol_weight * 1000
    except Exception:
        return 0.0


def calculate_esol(smiles: str) -> dict:
    """
    Calcula la solubilidad predicha usando el modelo ESOL de Delaney.
    
    Retorna:
        dict con logS_esol, solubility_class, mol_per_l, mg_per_l
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES inválido: {smiles}")

    # Descriptores necesarios para ESOL
    mw      = Descriptors.MolWt(mol)
    clogp   = Descriptors.MolLogP(mol)
    rb      = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    # Fracción de carbonos aromáticos (AP)
    total_c = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
    arom_c  = sum(1 for a in mol.GetAtoms() 
                  if a.GetAtomicNum() == 6 and a.GetIsAromatic())
    ap = arom_c / total_c if total_c > 0 else 0.0

    # Ecuación de Delaney
    logs_esol = 0.16 - 0.63*clogp - 0.0062*mw + 0.066*rb - 0.74*ap

    sol_class = classify_solubility(logs_esol)
    
    return {
        "logS_esol":         round(logs_esol, 3),
        "solubility_class":  sol_class["class"],
        "solubility_color":  sol_class["color"],
        "solubility_desc":   sol_class["description"],
        "mol_per_liter":     round(10 ** logs_esol, 6),
        "mg_per_liter":      round(logs_to_mg_per_l(logs_esol, mw), 2),
        "descriptors_used": {
            "MW":    round(mw, 2),
            "cLogP": round(clogp, 2),
            "RB":    rb,
            "AP":    round(ap, 3),
        }
    }
