"""
Flask REST API - RDKit Pharmacy Analysis Backend
Puerto: 5000
"""

import sys
import os
import base64
import json
import math
import io

# Agregar directorio actual al path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from flask import Flask, jsonify, request
from flask_cors import CORS
import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import (
    Draw, Descriptors, rdMolDescriptors, AllChem,
    DataStructs, Fragments
)
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import FilterCatalog
from PIL import Image

from esol_calculator import calculate_esol, classify_solubility

app = Flask(__name__)
CORS(app)

# ── Cargar dataset ESOL ─────────────────────────────────────────────────────
DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         '..', 'data', 'delaney_esol.csv')

_molecules_cache = None

def load_molecules():
    global _molecules_cache
    if _molecules_cache is not None:
        return _molecules_cache

    df = pd.read_csv(DATA_PATH)
    # Normalizar columnas
    col_map = {}
    for c in df.columns:
        cl = c.lower().strip()
        if 'smiles' in cl:
            col_map[c] = 'smiles'
        elif 'compound' in cl or 'name' in cl:
            col_map[c] = 'name'
        elif 'measured' in cl:
            col_map[c] = 'logs_measured'
        elif 'esol' in cl or 'predicted' in cl:
            col_map[c] = 'logs_esol'
    df = df.rename(columns=col_map)

    molecules = []
    for idx, row in df.iterrows():
        smiles = str(row.get('smiles', '')).strip()
        if not smiles or smiles == 'nan':
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        logs_m = float(row.get('logs_measured', 0)) if 'logs_measured' in row else None
        sol_class = classify_solubility(logs_m if logs_m is not None else -5)

        molecules.append({
            "id":             idx,
            "name":           str(row.get('name', f'Molecula_{idx}')).strip(),
            "smiles":         smiles,
            "logs_measured":  round(logs_m, 3) if logs_m is not None else None,
            "solubility_class":  sol_class["class"],
            "solubility_color":  sol_class["color"],
        })

    _molecules_cache = molecules
    return molecules


# ── Helpers ──────────────────────────────────────────────────────────────────

def mol_to_image_b64(mol, size=(400, 300)):
    """Genera imagen 2D de una molécula en base64."""
    drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    drawer.drawOptions().addStereoAnnotation = True
    drawer.drawOptions().addAtomIndices = False
    AllChem.Compute2DCoords(mol)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    bio = io.BytesIO(drawer.GetDrawingText())
    img = Image.open(bio)
    buf = io.BytesIO()
    img.save(buf, format='PNG')
    return base64.b64encode(buf.getvalue()).decode('utf-8')


def smiles_from_request():
    data = request.get_json(force=True)
    smiles = data.get('smiles', '').strip()
    if not smiles:
        raise ValueError("SMILES requerido")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES inválido: {smiles}")
    return smiles, mol


# ── Endpoints ────────────────────────────────────────────────────────────────

@app.route('/molecules', methods=['GET'])
def get_molecules():
    page     = int(request.args.get('page', 1))
    per_page = int(request.args.get('per_page', 50))
    search   = request.args.get('search', '').lower()
    sol_filter = request.args.get('solubility', '').strip()

    mols = load_molecules()

    if search:
        mols = [m for m in mols if search in m['name'].lower() or search in m['smiles'].lower()]
    if sol_filter:
        mols = [m for m in mols if m['solubility_class'] == sol_filter]

    total = len(mols)
    start = (page - 1) * per_page
    end   = start + per_page

    return jsonify({
        "total":     total,
        "page":      page,
        "per_page":  per_page,
        "molecules": mols[start:end]
    })


@app.route('/molecules/all', methods=['GET'])
def get_all_molecules():
    return jsonify(load_molecules())


@app.route('/molecule/<int:mol_id>', methods=['GET'])
def get_molecule(mol_id):
    mols = load_molecules()
    mol_data = next((m for m in mols if m['id'] == mol_id), None)
    if mol_data is None:
        return jsonify({"error": "Molécula no encontrada"}), 404
    return jsonify(mol_data)


@app.route('/analyze/structure', methods=['POST'])
def analyze_structure():
    try:
        smiles, mol = smiles_from_request()
        img_b64 = mol_to_image_b64(mol)
        
        inchi = Chem.MolToInchi(mol) or ''
        inchikey = Chem.InchiToInchiKey(inchi) if inchi else ''
        formula = rdMolDescriptors.CalcMolFormula(mol)
        
        return jsonify({
            "image_base64": img_b64,
            "smiles":       smiles,
            "inchi":        inchi,
            "inchikey":     inchikey,
            "formula":      formula,
        })
    except Exception as e:
        return jsonify({"error": str(e)}), 400


@app.route('/analyze/solubility', methods=['POST'])
def analyze_solubility():
    try:
        smiles, mol = smiles_from_request()
        data       = request.get_json(force=True)
        logs_meas  = data.get('logs_measured')

        esol = calculate_esol(smiles)

        result = {**esol}
        if logs_meas is not None:
            try:
                logs_m = float(logs_meas)
                sol_m  = classify_solubility(logs_m)
                result["logs_measured"]          = round(logs_m, 3)
                result["measured_class"]         = sol_m["class"]
                result["measured_color"]         = sol_m["color"]
                result["error_esol_vs_measured"] = round(abs(esol["logS_esol"] - logs_m), 3)
            except Exception:
                pass

        return jsonify(result)
    except Exception as e:
        return jsonify({"error": str(e)}), 400


@app.route('/analyze/properties', methods=['POST'])
def analyze_properties():
    try:
        smiles, mol = smiles_from_request()
        
        props = {
            "MW":              round(Descriptors.MolWt(mol), 2),
            "HeavyAtomMW":     round(Descriptors.HeavyAtomMolWt(mol), 2),
            "LogP":            round(Descriptors.MolLogP(mol), 3),
            "TPSA":            round(Descriptors.TPSA(mol), 2),
            "HBD":             rdMolDescriptors.CalcNumHBD(mol),
            "HBA":             rdMolDescriptors.CalcNumHBA(mol),
            "RotBonds":        rdMolDescriptors.CalcNumRotatableBonds(mol),
            "AromaticRings":   rdMolDescriptors.CalcNumAromaticRings(mol),
            "Rings":           rdMolDescriptors.CalcNumRings(mol),
            "HeavyAtoms":      mol.GetNumHeavyAtoms(),
            "Stereocenters":   len(Chem.FindMolChiralCenters(mol, includeUnassigned=True)),
            "FractionCSP3":    round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
            "MolarRefractivity": round(Descriptors.MolMR(mol), 2),
            "Complexity":      round(Descriptors.BertzCT(mol), 1),
        }
        
        # Contexto educativo para cada propiedad
        DESCRIPTIONS = {
            "MW":          ("Peso Molecular (Da)", "Masa total de la molécula. Lipinski: ≤ 500 Da"),
            "LogP":        ("LogP (Lipofilicidad)", "Coeficiente de partición octanol/agua. Lipinski: ≤ 5"),
            "TPSA":        ("Área Polar Superficial (Ų)", "Área ocupada por grupos polares. Absorción oral ↓ si > 140 Ų"),
            "HBD":         ("Donadores H-Bond", "Átomos N-H o O-H que donan puentes de H. Lipinski: ≤ 5"),
            "HBA":         ("Aceptores H-Bond", "Átomos N u O que aceptan puentes de H. Lipinski: ≤ 10"),
            "RotBonds":    ("Enlaces Rotables", "Flexibilidad molecular. > 10 puede afectar biodisponibilidad"),
            "AromaticRings":("Anillos Aromáticos", "Número de sistemas aromáticos. Relacionado con metabolismo"),
            "Rings":       ("Anillos Totales", "Ciclos en la estructura molecular"),
            "HeavyAtoms":  ("Átomos Pesados", "Número de átomos no hidrógeno"),
            "FractionCSP3":("Fracción Csp3", "Carbones sp3 / total. Relacionado con 3D-complexity"),
            "Complexity":  ("Complejidad Molecular", "Índice de Bertz-CT. Mide la complejidad topológica"),
        }
        
        result = []
        for key, value in props.items():
            label, desc = DESCRIPTIONS.get(key, (key, ""))
            result.append({"key": key, "label": label, "value": value, "description": desc})
        
        return jsonify({"properties": result})
    except Exception as e:
        return jsonify({"error": str(e)}), 400


@app.route('/analyze/lipinski', methods=['POST'])
def analyze_lipinski():
    try:
        smiles, mol = smiles_from_request()

        mw   = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd  = rdMolDescriptors.CalcNumHBD(mol)
        hba  = rdMolDescriptors.CalcNumHBA(mol)
        rb   = rdMolDescriptors.CalcNumRotatableBonds(mol)
        tpsa = Descriptors.TPSA(mol)

        rules = [
            {
                "rule": "Regla 1 — Peso Molecular",
                "criterion": "≤ 500 Da",
                "value": round(mw, 1),
                "unit": "Da",
                "pass": mw <= 500,
                "explanation": "Moléculas muy grandes (> 500 Da) tienen dificultad para atravesar membranas celulares."
            },
            {
                "rule": "Regla 2 — Lipofilicidad (LogP)",
                "criterion": "≤ 5",
                "value": round(logp, 2),
                "unit": "",
                "pass": logp <= 5,
                "explanation": "Un LogP muy alto indica poca solubilidad en agua e incapacidad de disolverse en fluidos corporales."
            },
            {
                "rule": "Regla 3 — Donadores H-Bond",
                "criterion": "≤ 5",
                "value": hbd,
                "unit": "",
                "pass": hbd <= 5,
                "explanation": "Demasiados donadores de hidrógeno dificultan la absorción a través de membranas lipídicas."
            },
            {
                "rule": "Regla 4 — Aceptores H-Bond",
                "criterion": "≤ 10",
                "value": hba,
                "unit": "",
                "pass": hba <= 10,
                "explanation": "Un exceso de aceptores reduce la permeabilidad de la membrana intestinal."
            },
        ]

        # Reglas adicionales de Veber (2002)
        veber_rules = [
            {
                "rule": "Veber — TPSA",
                "criterion": "≤ 140 Ų",
                "value": round(tpsa, 1),
                "unit": "Ų",
                "pass": tpsa <= 140,
                "explanation": "Un TPSA ≤ 140 Ų predice buena permeabilidad oral (Veber, 2002)."
            },
            {
                "rule": "Veber — enlaces Rotables",
                "criterion": "≤ 10",
                "value": rb,
                "unit": "",
                "pass": rb <= 10,
                "explanation": "Pocas rotaciones implica conformaciones más rígidas y mejor absorción."
            },
        ]

        passes = sum(1 for r in rules if r["pass"])
        failures = 4 - passes
        drug_like = failures <= 1  # Lipinski permite 1 violación

        return jsonify({
            "lipinski_rules": rules,
            "veber_rules":    veber_rules,
            "passes":         passes,
            "failures":       failures,
            "drug_like":      drug_like,
            "verdict":        "✅ Drug-like (cumple Lipinski)" if drug_like else "❌ No drug-like (viola ≥ 2 reglas)",
        })
    except Exception as e:
        return jsonify({"error": str(e)}), 400


@app.route('/analyze/admet', methods=['POST'])
def analyze_admet():
    try:
        smiles, mol = smiles_from_request()

        mw   = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd  = rdMolDescriptors.CalcNumHBD(mol)
        hba  = rdMolDescriptors.CalcNumHBA(mol)
        tpsa = Descriptors.TPSA(mol)
        rb   = rdMolDescriptors.CalcNumRotatableBonds(mol)

        # Puntuaciones ADMET (0-100, heurísticas basadas en propiedades)
        def score_absorption():
            # Basado en TPSA y HBD (regla de Veber)
            s = 100
            if tpsa > 140: s -= 50
            elif tpsa > 90: s -= 20
            if hbd > 5: s -= 30
            if mw > 500: s -= 20
            return max(0, min(100, s))

        def score_distribution():
            # LogP óptimo 1-3 para BBB; HBA baja mejora distribución
            s = 70
            if 1 <= logp <= 3: s += 20
            elif logp < 0 or logp > 5: s -= 30
            if hba > 10: s -= 20
            return max(0, min(100, s))

        def score_metabolism():
            # Pocos anillos aromáticos y grupos funcionales reactivos = mejor metabolismo
            ar = rdMolDescriptors.CalcNumAromaticRings(mol)
            s = 80
            if ar > 3: s -= 25
            if logp > 4: s -= 15
            return max(0, min(100, s))

        def score_excretion():
            # MW y TPSA influyen en excreción biliar vs renal
            s = 75
            if mw < 300: s += 15
            elif mw > 500: s -= 20
            if tpsa > 100: s -= 10
            return max(0, min(100, s))

        def score_toxicity():
            # Heurístico inverso: menos anillos aromáticos, menos grupos reactivos = menos tóxico
            ar = rdMolDescriptors.CalcNumAromaticRings(mol)
            s = 80
            if ar > 2: s -= ar * 5
            if logp > 4.5: s -= 15
            num_nitro = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 7 and a.GetNoImplicit())
            if num_nitro > 2: s -= 20
            return max(0, min(100, s))

        scores = {
            "Absorción":     score_absorption(),
            "Distribución":  score_distribution(),
            "Metabolismo":   score_metabolism(),
            "Excreción":     score_excretion(),
            "Toxicidad":     score_toxicity(),  # invertido: 100 = no tóxico
        }

        descriptions = {
            "Absorción":    "Capacidad de la molécula para ser absorbida en el tracto GI",
            "Distribución": "Distribución en tejidos y potencial para cruzar la barrera hematoencefálica",
            "Metabolismo":  "Susceptibilidad al metabolismo hepático (CYP450)",
            "Excreción":    "Velocidad y vía de eliminación del organismo",
            "Toxicidad":    "Perfil de seguridad estimado (100 = muy segura, 0 = alta toxicidad)",
        }

        admet_list = [
            {"axis": k, "score": v, "description": descriptions[k]}
            for k, v in scores.items()
        ]

        overall = round(sum(scores.values()) / len(scores), 1)

        return jsonify({
            "admet": admet_list,
            "overall_score": overall,
            "interpretation": (
                "Excelente perfil ADMET" if overall >= 75 else
                "Buen perfil ADMET" if overall >= 60 else
                "Perfil ADMET moderado" if overall >= 45 else
                "Perfil ADMET bajo — considerar optimización"
            )
        })
    except Exception as e:
        return jsonify({"error": str(e)}), 400


@app.route('/analyze/similarity', methods=['POST'])
def analyze_similarity():
    try:
        smiles, mol = smiles_from_request()
        data    = request.get_json(force=True)
        top_n   = int(data.get('top_n', 10))
        radius  = int(data.get('radius', 2))
        nbits   = int(data.get('nbits', 2048))

        query_fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)

        all_mols  = load_molecules()
        results   = []

        for m in all_mols:
            if m['smiles'] == smiles:
                continue
            mol2 = Chem.MolFromSmiles(m['smiles'])
            if mol2 is None:
                continue
            fp2  = AllChem.GetMorganFingerprintAsBitVect(mol2, radius, nBits=nbits)
            sim  = DataStructs.TanimotoSimilarity(query_fp, fp2)
            results.append({
                "id":         m["id"],
                "name":       m["name"],
                "smiles":     m["smiles"],
                "similarity": round(sim, 4),
                "logs_measured": m.get("logs_measured"),
                "solubility_class": m.get("solubility_class"),
                "solubility_color": m.get("solubility_color"),
            })

        results.sort(key=lambda x: x["similarity"], reverse=True)
        return jsonify({"similar_molecules": results[:top_n]})
    except Exception as e:
        return jsonify({"error": str(e)}), 400


@app.route('/analyze/descriptors', methods=['POST'])
def analyze_descriptors():
    try:
        smiles, mol = smiles_from_request()

        desc_list = []
        for name, fn in Descriptors.descList:
            try:
                val = fn(mol)
                if isinstance(val, float):
                    if math.isnan(val) or math.isinf(val):
                        val = None
                    else:
                        val = round(val, 4)
                desc_list.append({"name": name, "value": val})
            except Exception:
                desc_list.append({"name": name, "value": None})

        return jsonify({
            "total": len(desc_list),
            "descriptors": desc_list
        })
    except Exception as e:
        return jsonify({"error": str(e)}), 400


@app.route('/analyze/substructure', methods=['POST'])
def analyze_substructure():
    try:
        data   = request.get_json(force=True)
        smiles = data.get('smiles', '').strip()
        smarts = data.get('smarts', '').strip()

        if not smarts:
            return jsonify({"error": "SMARTS requerido"}), 400

        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            return jsonify({"error": "SMARTS inválido"}), 400

        all_mols = load_molecules()
        matches  = []
        for m in all_mols:
            mol2 = Chem.MolFromSmiles(m['smiles'])
            if mol2 and mol2.HasSubstructMatch(pattern):
                matches.append({
                    "id":    m["id"],
                    "name":  m["name"],
                    "smiles": m["smiles"],
                    "solubility_class": m.get("solubility_class"),
                    "solubility_color": m.get("solubility_color"),
                })

        return jsonify({
            "smarts":  smarts,
            "matches": matches,
            "count":   len(matches)
        })
    except Exception as e:
        return jsonify({"error": str(e)}), 400


@app.route('/molecule/validate', methods=['POST'])
def validate_molecule():
    """Valida un SMILES y devuelve imagen + propiedades básicas + Lipinski."""
    try:
        data   = request.get_json(force=True)
        smiles = data.get('smiles', '').strip()
        name   = data.get('name', 'Sin nombre').strip() or 'Sin nombre'

        if not smiles:
            return jsonify({"valid": False, "error": "SMILES vacío"}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({
                "valid":   False,
                "error":   "SMILES inválido: RDKit no puede interpretar esta cadena.",
                "verdict": "❌ No es una molécula válida"
            }), 200

        # Imagen 2D
        img_b64 = mol_to_image_b64(mol, size=(350, 260))

        # Propiedades básicas
        mw   = round(Descriptors.MolWt(mol), 2)
        logp = round(Descriptors.MolLogP(mol), 3)
        hbd  = rdMolDescriptors.CalcNumHBD(mol)
        hba  = rdMolDescriptors.CalcNumHBA(mol)
        tpsa = round(Descriptors.TPSA(mol), 2)
        rb   = rdMolDescriptors.CalcNumRotatableBonds(mol)
        formula = rdMolDescriptors.CalcMolFormula(mol)

        # Lipinski simplificado
        violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
        drug_like  = violations <= 1

        # ESOL
        from esol_calculator import calculate_esol, classify_solubility
        esol = calculate_esol(smiles)
        sol_class = esol.get('solubility_class', '—')
        sol_color = esol.get('solubility_color', '#94a3b8')
        logs_pred = esol.get('logS_esol', None)

        return jsonify({
            "valid":       True,
            "name":        name,
            "smiles":      Chem.MolToSmiles(mol),
            "formula":     formula,
            "image_base64": img_b64,
            "properties": {
                "MW":      mw,
                "LogP":    logp,
                "HBD":     hbd,
                "HBA":     hba,
                "TPSA":    tpsa,
                "RotBonds": rb,
            },
            "lipinski": {
                "violations": violations,
                "drug_like":  drug_like,
                "verdict":    "✅ Drug-like (Lipinski OK)" if drug_like else f"⚠️ No drug-like ({violations} violaciones)",
            },
            "solubility": {
                "logS_esol":       logs_pred,
                "solubility_class": sol_class,
                "solubility_color": sol_color,
            }
        })
    except Exception as e:
        return jsonify({"valid": False, "error": str(e)}), 500


@app.route('/molecule/add', methods=['POST'])
def add_molecule():
    """Agrega una molécula al CSV y recarga la caché."""
    global _molecules_cache
    try:
        data   = request.get_json(force=True)
        smiles = data.get('smiles', '').strip()
        name   = data.get('name', '').strip()
        logs_m = data.get('logs_measured', None)

        if not smiles or not name:
            return jsonify({"error": "SMILES y nombre son requeridos"}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({"error": "SMILES inválido"}), 400

        canonical = Chem.MolToSmiles(mol)

        # Verificar duplicado por SMILES canónico
        existing = load_molecules()
        for m in existing:
            if Chem.MolToSmiles(Chem.MolFromSmiles(m['smiles'])) == canonical:
                return jsonify({"error": f"La molécula ya existe en la biblioteca como '{m['name']}'"}), 409

        # Calcular ESOL
        from esol_calculator import calculate_esol
        esol = calculate_esol(canonical)
        logs_pred = esol.get('logS_esol', 0)

        if logs_m is not None:
            try:
                logs_m = float(logs_m)
            except Exception:
                logs_m = logs_pred
        else:
            logs_m = logs_pred

        # Escribir al CSV
        import csv as csv_mod
        data_path = os.path.normpath(DATA_PATH)
        with open(data_path, 'a', newline='', encoding='utf-8') as f:
            writer = csv_mod.writer(f)
            writer.writerow([name, canonical, round(logs_m, 3), round(logs_pred, 3)])

        # Invalidar caché para recargar
        _molecules_cache = None

        # Devolver la nueva molécula con su ID
        mols = load_molecules()
        new_mol = next((m for m in mols if m['smiles'] == canonical), None)

        return jsonify({
            "success": True,
            "message": f"'{name}' agregada a la biblioteca.",
            "molecule": new_mol
        })
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route('/molecule/delete', methods=['DELETE'])
def delete_molecule():
    """Elimina permanentemente una molécula del CSV usando su SMILES canónico."""
    global _molecules_cache
    try:
        data   = request.get_json(force=True)
        smiles = data.get('smiles', '').strip()

        if not smiles:
            return jsonify({"error": "SMILES es requerido para eliminar"}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({"error": "SMILES inválido"}), 400

        canonical = Chem.MolToSmiles(mol)
        
        # Leer el CSV actual
        import csv as csv_mod
        data_path = os.path.normpath(DATA_PATH)
        
        rows_to_keep = []
        deleted = False
        
        with open(data_path, 'r', encoding='utf-8') as f:
            reader = csv_mod.reader(f)
            header = next(reader, None)
            if header:
                rows_to_keep.append(header)
            
            for row in reader:
                if len(row) >= 2:
                    row_smiles = row[1]
                    try:
                        row_mol = Chem.MolFromSmiles(row_smiles)
                        row_canonical = Chem.MolToSmiles(row_mol) if row_mol else row_smiles
                    except Exception:
                        row_canonical = row_smiles
                        
                    if row_canonical == canonical:
                        deleted = True
                        continue  # Skipear esta fila
                rows_to_keep.append(row)

        if not deleted:
            return jsonify({"error": "Molécula no encontrada en la biblioteca."}), 404

        # Sobreescribir el CSV con las filas restantes
        with open(data_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv_mod.writer(f)
            writer.writerows(rows_to_keep)

        # Invalidar caché para recargar en la próxima petición
        _molecules_cache = None

        return jsonify({
            "success": True,
            "message": "Molécula eliminada permanentemente"
        })
    except Exception as e:
        return jsonify({"error": str(e)}), 500


# ── Reacciones predefinidas (SMARTS) ─────────────────────────────────────────

REACTIONS = {
    "esterification": {
        "name": "Esterificación",
        "smarts": "[C:1](=O)[OH:2].[O:3][C:4]>>[C:1](=O)[O:3][C:4]",
        "needs_two": True,
        "r1_hint": "Ácido carboxílico (ej: CC(=O)O)",
        "r2_hint": "Alcohol (ej: CCO)",
        "explanation": "Un ácido carboxílico reacciona con un alcohol formando un éster y eliminando agua. Reacción clave en síntesis de fragancias y medicamentos.",
    },
    "amide": {
        "name": "Formación de Amida",
        "smarts": "[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[N:2]",
        "needs_two": True,
        "r1_hint": "Ácido carboxílico (ej: CC(=O)O)",
        "r2_hint": "Amina primaria (ej: CN)",
        "explanation": "Un ácido carboxílico reacciona con una amina primaria para formar un enlace amida (péptido). Fundamental en bioquímica y síntesis de fármacos.",
    },
    "imine": {
        "name": "Formación de Imina (Condensación)",
        "smarts": "[C:1]=O.[NH2:2]>>[C:1]=[N:2]",
        "needs_two": True,
        "r1_hint": "Aldehído o cetona (ej: CC=O)",
        "r2_hint": "Amina primaria (ej: N)",
        "explanation": "Un carbonilo (aldehído o cetona) reacciona con una amina primaria formando una imina (base de Schiff) con eliminación de agua.",
    },
    "reduction_ketone": {
        "name": "Reducción de Cetona/Aldehído",
        "smarts": "[C:1](=[O:2])[C:3]>>[C:1]([OH])[C:3]",
        "needs_two": False,
        "r1_hint": "Cetona (ej: CC(=O)C) o aldehído (ej: CC=O)",
        "r2_hint": "",
        "explanation": "La reducción de un carbonilo produce un alcohol. Simula la adición de H₂ (hidrogenación catalítica) o uso de NaBH₄/LiAlH₄.",
    },
    "oxidation_alcohol": {
        "name": "Oxidación de Alcohol Primario",
        "smarts": "[C:1][CH2:2][OH]>>[C:1][C:2]=O",
        "needs_two": False,
        "r1_hint": "Alcohol primario (ej: CCO, CCCO)",
        "r2_hint": "",
        "explanation": "Los alcoholes primarios se oxidan a aldehídos con agentes como PCC, DMSO activado (Swern) o MnO₂. Reacción de uso frecuente en síntesis.",
    },
    "alkylation": {
        "name": "N-Alquilación (Sustitución)",
        "smarts": "[NH2:1].[CH2:2][Br]>>[N:1][C:2]",
        "needs_two": True,
        "r1_hint": "Amina primaria (ej: CN, CCN)",
        "r2_hint": "Haluro de alquilo (ej: CBr)",
        "explanation": "Una amina primaria reacciona con un haluro de alquilo (SN2) formando una amina secundaria. Base sintética de muchos fármacos y colorantes.",
    },
    "halogenation": {
        "name": "Halogenación Aromática (EAS)",
        "smarts": "[c:1][H]>>[c:1]Cl",
        "needs_two": False,
        "r1_hint": "Aromático (ej: c1ccccc1, Cc1ccccc1)",
        "r2_hint": "",
        "explanation": "Sustitución Electrófila Aromática (SEA). El cloro reemplaza un H del anillo aromático en presencia de un catalizador Lewis (AlCl₃ o FeCl₃).",
    },
    "diels_alder": {
        "name": "Cicloadición Diels-Alder [4+2]",
        "smarts": "[C:1]=[C:2][C:3]=[C:4].[C:5]=[C:6]>>[C:1]1[C:2]=[C:3][C:4][C:6][C:5]1",
        "needs_two": True,
        "r1_hint": "Dieno (ej: C=CC=C, butadieno)",
        "r2_hint": "Dienófilo (ej: C=C, C=CC=O)",
        "explanation": "El dieno (4π) reacciona con el dienófilo (2π) en una cicloadición pericíclica [4+2] para formar un ciclo de 6 miembros. Reacción estereoespecífica.",
    },
}


@app.route('/analyze/reaction', methods=['POST'])
def analyze_reaction():
    """Simula una reacción química y devuelve los productos con RDKit."""
    try:
        from rdkit.Chem import rdChemReactions

        data     = request.get_json(force=True)
        rxn_key  = data.get('reaction_type', '').strip()
        smiles1  = data.get('reactant1', '').strip()
        smiles2  = data.get('reactant2', '').strip()

        if rxn_key not in REACTIONS:
            return jsonify({"error": f"Tipo de reacción desconocido: {rxn_key}"}), 400

        rxn_info = REACTIONS[rxn_key]

        if not smiles1:
            return jsonify({"error": "Se requiere al menos el Reactivo 1 (SMILES)"}), 400

        mol1 = Chem.MolFromSmiles(smiles1)
        if mol1 is None:
            return jsonify({"error": f"SMILES inválido para Reactivo 1: {smiles1}"}), 400

        mol2 = None
        if rxn_info["needs_two"]:
            if not smiles2:
                return jsonify({"error": "Esta reacción requiere dos reactivos (Reactivo 2 vacío)"}), 400
            mol2 = Chem.MolFromSmiles(smiles2)
            if mol2 is None:
                return jsonify({"error": f"SMILES inválido para Reactivo 2: {smiles2}"}), 400

        # Aplicar la reacción
        rxn = rdChemReactions.ReactionFromSmarts(rxn_info["smarts"])
        if rxn is None:
            return jsonify({"error": "Error interno: plantilla SMARTS inválida"}), 500

        reactants = (mol1, mol2) if mol2 else (mol1,)
        try:
            products_sets = rxn.RunReactants(reactants)
        except Exception as e:
            return jsonify({"error": f"La reacción no pudo ejecutarse: {str(e)}"}), 400

        if not products_sets:
            return jsonify({
                "error": "No se generaron productos. Verifica que los reactivos coincidan con el tipo de reacción seleccionado.",
                "hint": f"R1 esperado: {rxn_info['r1_hint']}" + (f" | R2 esperado: {rxn_info['r2_hint']}" if rxn_info['needs_two'] else ""),
            }), 200

        # Recopilar productos únicos (SMILES canónico)
        seen = set()
        products_out = []
        for pset in products_sets:
            for pmol in pset:
                try:
                    Chem.SanitizeMol(pmol)
                    canon = Chem.MolToSmiles(pmol)
                    if canon in seen:
                        continue
                    seen.add(canon)

                    img_b64 = mol_to_image_b64(pmol, size=(320, 240))
                    mw   = round(Descriptors.MolWt(pmol), 2)
                    logp = round(Descriptors.MolLogP(pmol), 3)
                    tpsa = round(Descriptors.TPSA(pmol), 2)
                    hbd  = rdMolDescriptors.CalcNumHBD(pmol)
                    hba  = rdMolDescriptors.CalcNumHBA(pmol)
                    formula = rdMolDescriptors.CalcMolFormula(pmol)
                    violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])

                    from esol_calculator import calculate_esol
                    try:
                        esol = calculate_esol(canon)
                        logs_pred = esol.get('logS_esol')
                        sol_class = esol.get('solubility_class', '—')
                        sol_color = esol.get('solubility_color', '#94a3b8')
                    except Exception:
                        logs_pred = None
                        sol_class = '—'
                        sol_color = '#94a3b8'

                    products_out.append({
                        "smiles":   canon,
                        "formula":  formula,
                        "image_base64": img_b64,
                        "properties": {
                            "MW": mw, "LogP": logp, "TPSA": tpsa,
                            "HBD": hbd, "HBA": hba,
                        },
                        "lipinski": {
                            "violations": violations,
                            "drug_like": violations <= 1,
                            "verdict": "✅ Drug-like" if violations <= 1 else f"⚠️ {violations} violaciones",
                        },
                        "solubility": {
                            "logS_esol": logs_pred,
                            "class": sol_class,
                            "color": sol_color,
                        },
                    })
                except Exception:
                    continue

            if len(products_out) >= 5:
                break  # Limitar a 5 productos únicos

        if not products_out:
            return jsonify({"error": "Los productos generados no pudieron sanitizarse. Prueba otros reactivos."}), 200

        # Imágenes de los reactivos
        def mol_img(m):
            try:
                return mol_to_image_b64(m, size=(280, 200))
            except Exception:
                return None

        return jsonify({
            "reaction_name":  rxn_info["name"],
            "reaction_smarts": rxn_info["smarts"],
            "explanation":    rxn_info["explanation"],
            "needs_two":      rxn_info["needs_two"],
            "reactant1": {
                "smiles": Chem.MolToSmiles(mol1),
                "formula": rdMolDescriptors.CalcMolFormula(mol1),
                "image_base64": mol_img(mol1),
            },
            "reactant2": {
                "smiles": Chem.MolToSmiles(mol2) if mol2 else None,
                "formula": rdMolDescriptors.CalcMolFormula(mol2) if mol2 else None,
                "image_base64": mol_img(mol2) if mol2 else None,
            } if mol2 else None,
            "products": products_out,
        })
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route('/reactions', methods=['GET'])
def get_reactions():
    """Devuelve la lista de reacciones predefinidas disponibles."""
    return jsonify([
        {
            "key": k,
            "name": v["name"],
            "needs_two": v["needs_two"],
            "r1_hint": v["r1_hint"],
            "r2_hint": v["r2_hint"],
            "explanation": v["explanation"],
        }
        for k, v in REACTIONS.items()
    ])


@app.route('/health', methods=['GET'])
def health():
    return jsonify({"status": "ok", "rdkit": "loaded"})

@app.route('/shutdown', methods=['POST'])
def shutdown():
    """Apaga el servidor Flask de forma limpia y abrupta."""
    import os
    import threading
    def kill_server():
        os._exit(0)
    # Esperamos 300ms para asegurar que responda el HTTP 200 al frontend
    threading.Timer(0.3, kill_server).start()
    return jsonify({"success": True, "message": "Backend apagándose"})


if __name__ == '__main__':
    print("🚀 Iniciando servidor RDKit en http://localhost:5000")
    print("📂 Dataset ESOL cargando desde:", DATA_PATH)
    app.run(host='0.0.0.0', port=5000, debug=False)
