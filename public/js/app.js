/**
 * RDKit Pharma Lab — Frontend JavaScript
 * Maneja la interfaz, llamadas API y visualización de resultados
 */

const API = '/api';
let currentMolecule = null;
let allMolecules     = [];
let radarChart       = null;
let allDescriptors   = [];

// ─── Inicialización ───────────────────────────────────────────────────────────
document.addEventListener('DOMContentLoaded', async () => {
    setupTabs();
    setupFilters();
    setupDescSearch();
    setupSubstructure();
    setupMoleculeBuilder();
    setupReactionSimulator();
    await loadMoleculeLibrary();
    await loadReactionMeta();
});

// ─── Cargar biblioteca de moléculas ──────────────────────────────────────────
async function loadMoleculeLibrary() {
    try {
        const res  = await fetch(`${API}/molecules/all`);
        const data = await res.json();
        allMolecules = data;

        renderMolList(allMolecules);
        updateGlobalStats(allMolecules);

        document.getElementById('dbStatus').textContent =
            `✅ ${allMolecules.length} compuestos cargados`;
    } catch (err) {
        document.getElementById('dbStatus').textContent = '❌ Error al conectar con Python';
        document.getElementById('molList').innerHTML =
            `<div class="error-box" style="margin:12px;">
              ⚠️ No se puede conectar con el servidor Python.<br>
              <small>Asegúrate de que <strong>python/server.py</strong> esté corriendo.</small>
             </div>`;
    }
}

function updateGlobalStats(mols) {
    const total    = mols.length;
    const soluble  = mols.filter(m => ['Soluble','Muy soluble','Altamente soluble','Moderadamente soluble'].includes(m.solubility_class)).length;
    const insol    = mols.filter(m => m.solubility_class === 'Insoluble').length;
    const avgLogs  = mols.reduce((s,m) => s + (m.logs_measured ?? 0), 0) / total;

    document.getElementById('statTotal').textContent    = total;
    document.getElementById('statSoluble').textContent  = soluble;
    document.getElementById('statInsoluble').textContent = insol;
    document.getElementById('statAvgLogs').textContent  = avgLogs.toFixed(2);
}

// ─── Renderizar lista de moléculas ───────────────────────────────────────────
function renderMolList(mols) {
    const list = document.getElementById('molList');
    document.getElementById('molCount').textContent = `${mols.length} moléculas`;

    if (mols.length === 0) {
        list.innerHTML = '<div style="padding:20px;text-align:center;color:var(--text-muted);font-size:0.82rem;">Sin resultados</div>';
        return;
    }

    list.innerHTML = mols.map(m => `
        <div class="mol-item ${currentMolecule?.id === m.id ? 'selected' : ''}"
             data-id="${m.id}" onclick="selectMolecule(${m.id})">
          <div class="sol-dot" style="background:${m.solubility_color}"></div>
          <div class="mol-name" title="${m.name}">${m.name}</div>
          <div class="mol-logs">${m.logs_measured !== null ? m.logs_measured.toFixed(2) : '—'}</div>
        </div>
    `).join('');
}

// ─── Filtros y búsqueda ──────────────────────────────────────────────────────
function setupFilters() {
    let activeFilter = '';
    document.querySelectorAll('.filter-btn').forEach(btn => {
        btn.addEventListener('click', () => {
            document.querySelectorAll('.filter-btn').forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
            activeFilter = btn.dataset.filter;
            applyFilters();
        });
    });

    document.getElementById('molSearch').addEventListener('input', applyFilters);

    function applyFilters() {
        const q = document.getElementById('molSearch').value.toLowerCase();
        let mols = allMolecules;
        if (activeFilter) mols = mols.filter(m => m.solubility_class === activeFilter);
        if (q) mols = mols.filter(m => m.name.toLowerCase().includes(q) || m.smiles.toLowerCase().includes(q));
        renderMolList(mols);
    }
}

// ─── Seleccionar molécula ────────────────────────────────────────────────────
async function selectMolecule(id) {
    const mol = allMolecules.find(m => m.id === id);
    if (!mol) return;

    currentMolecule = mol;

    // Resaltar en sidebar
    document.querySelectorAll('.mol-item').forEach(el => {
        el.classList.toggle('selected', parseInt(el.dataset.id) === id);
    });

    // Mostrar panel de análisis
    document.getElementById('welcomeScreen').style.display  = 'none';
    document.getElementById('analysisPanel').style.display  = 'block';

    // Actualizar Header info
    document.getElementById('panelMolName').textContent = mol.name;
    document.getElementById('panelSmiles').textContent  = mol.smiles;
    
    // Show delete button
    const deleteBtn = document.getElementById('btnDeleteMol');
    if (deleteBtn) deleteBtn.style.display = 'flex';
    const badge = document.getElementById('panelSolBadge');
    badge.textContent = mol.solubility_class;
    badge.style.background = hexToRgba(mol.solubility_color, 0.2);
    badge.style.border      = `1px solid ${hexToRgba(mol.solubility_color, 0.5)}`;
    badge.style.color       = mol.solubility_color;

    // Cargar pestaña activa
    const activeTab = document.querySelector('.tab-btn.active')?.dataset.tab || 'estructura';
    loadTab(activeTab, mol);
}

// ─── Tabs ────────────────────────────────────────────────────────────────────
function setupTabs() {
    document.querySelectorAll('.tab-btn').forEach(btn => {
        btn.addEventListener('click', () => {
            document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
            document.querySelectorAll('.tab-panel').forEach(p => p.classList.remove('active'));
            btn.classList.add('active');
            document.getElementById(`tab-${btn.dataset.tab}`).classList.add('active');
            if (currentMolecule) loadTab(btn.dataset.tab, currentMolecule);
        });
    });
}

function loadTab(tab, mol) {
    switch (tab) {
        case 'estructura':   loadStructure(mol);    break;
        case 'solubilidad':  loadSolubility(mol);   break;
        case 'propiedades':  loadProperties(mol);   break;
        case 'lipinski':     loadLipinski(mol);     break;
        case 'admet':        loadADMET(mol);         break;
        case 'similitud':    loadSimilarity(mol);   break;
        case 'descriptores': loadDescriptors(mol);  break;
        case 'subestructura': break; // manual search
    }
}

// ─── Tab: Estructura ─────────────────────────────────────────────────────────
async function loadStructure(mol) {
    const el = document.getElementById('structureContent');
    el.innerHTML = spinner();
    try {
        const res  = await postAPI('/analyze/structure', { smiles: mol.smiles });
        const data = await res.json();
        el.innerHTML = `
            <div class="mol-image-wrap">
              <img src="data:image/png;base64,${data.image_base64}" alt="Estructura 2D de ${mol.name}">
            </div>
            <div class="info-grid">
              <div class="info-card">
                <div class="info-card-label">Fórmula</div>
                <div class="info-card-value" style="font-size:1rem;">${data.formula}</div>
              </div>
              <div class="info-card" style="grid-column:1/-1;">
                <div class="info-card-label">SMILES</div>
                <div style="font-family:'Courier New',monospace;font-size:0.78rem;color:var(--accent-teal);word-break:break-all;">${data.smiles}</div>
              </div>
              <div class="info-card" style="grid-column:1/-1;">
                <div class="info-card-label">InChI</div>
                <div style="font-family:'Courier New',monospace;font-size:0.72rem;color:var(--text-muted);word-break:break-all;">${data.inchi || '—'}</div>
              </div>
              <div class="info-card" style="grid-column:1/-1;">
                <div class="info-card-label">InChIKey</div>
                <div style="font-family:'Courier New',monospace;font-size:0.78rem;color:var(--accent-purple);">${data.inchikey || '—'}</div>
              </div>
            </div>`;
    } catch (e) {
        el.innerHTML = errorBox(e);
    }
}

// ─── Tab: Solubilidad ─────────────────────────────────────────────────────────
async function loadSolubility(mol) {
    const el = document.getElementById('solubilityContent');
    el.innerHTML = spinner();
    try {
        const res  = await postAPI('/analyze/solubility', {
            smiles: mol.smiles,
            logs_measured: mol.logs_measured
        });
        const d = await res.json();

        const logsRange = { min: -12, max: 2 };
        const pct = v => Math.max(0, Math.min(100,
            ((v - logsRange.min) / (logsRange.max - logsRange.min)) * 100));

        const measuredRow = d.logs_measured !== undefined ? `
            <div class="comparison-row">
                <div class="cbar-label">Medido</div>
                <div class="cbar-track"><div class="cbar-fill" style="width:${pct(d.logs_measured)}%;background:linear-gradient(90deg,#68d391,#38a169)"></div></div>
                <div class="cbar-val" style="color:#68d391">${d.logs_measured}</div>
            </div>
            <div class="comparison-row">
                <div class="cbar-label">Error |Δ|</div>
                <div class="cbar-track"><div class="cbar-fill" style="width:${Math.min(50,d.error_esol_vs_measured)*10}%;background:linear-gradient(90deg,#f6ad55,#e53e3e)"></div></div>
                <div class="cbar-val" style="color:#f6ad55">${d.error_esol_vs_measured}</div>
            </div>` : '';

        el.innerHTML = `
            <div style="display:grid;grid-template-columns:1fr 1fr;gap:16px;align-items:start;">
              <div class="logs-display">
                <div class="logs-label">logS (ESOL predicho)</div>
                <div class="logs-value">${d.logS_esol}</div>
                <div class="logs-label">log(mol/L)</div>
                <div class="logs-class" style="background:${hexToRgba(d.solubility_color,0.2)};border:1px solid ${hexToRgba(d.solubility_color,0.5)};color:${d.solubility_color};">${d.solubility_class}</div>
                <div style="font-size:0.75rem;color:var(--text-muted);margin-top:8px;">${d.solubility_desc}</div>
              </div>
              <div>
                <div style="font-size:0.75rem;color:var(--text-muted);margin-bottom:8px;">Conversiones</div>
                <div class="info-card" style="margin-bottom:8px;">
                  <div class="info-card-label">mol/L</div>
                  <div class="info-card-value" style="font-size:1rem;">${d.mol_per_liter.toExponential(2)}</div>
                </div>
                <div class="info-card">
                  <div class="info-card-label">mg/L</div>
                  <div class="info-card-value" style="font-size:1rem;">${d.mg_per_liter.toLocaleString('es')}</div>
                </div>
              </div>
            </div>
            <div style="margin-top:16px;">
              <div style="font-size:0.75rem;color:var(--text-muted);margin-bottom:10px;text-transform:uppercase;letter-spacing:.05em;">Comparación logS</div>
              <div class="comparison-row">
                <div class="cbar-label">ESOL pred.</div>
                <div class="cbar-track"><div class="cbar-fill" style="width:${pct(d.logS_esol)}%;background:linear-gradient(90deg,#4facfe,#00f2fe)"></div></div>
                <div class="cbar-val" style="color:var(--accent-blue)">${d.logS_esol}</div>
              </div>
              ${measuredRow}
            </div>
            <div style="margin-top:16px;border-top:1px solid var(--border);padding-top:14px;">
              <div style="font-size:0.75rem;color:var(--text-muted);margin-bottom:8px;">Descriptores usados en ESOL</div>
              <div class="info-grid">
                ${Object.entries(d.descriptors_used).map(([k,v]) =>
                    `<div class="info-card"><div class="info-card-label">${k}</div><div class="info-card-value" style="font-size:1rem;">${v}</div></div>`
                ).join('')}
              </div>
            </div>`;
    } catch (e) {
        el.innerHTML = errorBox(e);
    }
}

// ─── Tab: Propiedades ────────────────────────────────────────────────────────
async function loadProperties(mol) {
    const el = document.getElementById('propertiesContent');
    el.innerHTML = spinner();
    try {
        const res  = await postAPI('/analyze/properties', { smiles: mol.smiles });
        const data = await res.json();
        el.innerHTML = `<div class="info-grid">
            ${data.properties.map(p => `
                <div class="info-card">
                  <div class="info-card-label">${p.label}</div>
                  <div class="info-card-value">${typeof p.value === 'number' ? p.value : '—'}</div>
                  <div class="info-card-desc">${p.description}</div>
                </div>`).join('')}
            </div>`;
    } catch (e) {
        el.innerHTML = errorBox(e);
    }
}

// ─── Tab: Lipinski ───────────────────────────────────────────────────────────
async function loadLipinski(mol) {
    const el = document.getElementById('lipinskiContent');
    el.innerHTML = spinner();
    try {
        const res  = await postAPI('/analyze/lipinski', { smiles: mol.smiles });
        const d    = await res.json();

        const ruleHtml = r => `
            <div class="rule-row ${r.pass ? 'pass' : 'fail'}">
              <div class="rule-icon">${r.pass ? '✅' : '❌'}</div>
              <div class="rule-info">
                <div class="rule-name">${r.rule}</div>
                <div class="rule-explanation">${r.explanation}</div>
              </div>
              <div class="rule-value">
                <div class="rule-val-num" style="color:${r.pass ? '#68d391' : '#fc8181'}">${r.value}${r.unit ? ' ' + r.unit : ''}</div>
                <div class="rule-criterion">${r.criterion}</div>
              </div>
            </div>`;

        el.innerHTML = `
            <div class="lipinski-grid">
              <div style="font-size:0.78rem;color:var(--text-muted);margin-bottom:6px;font-weight:600;">REGLAS DE LIPINSKI</div>
              ${d.lipinski_rules.map(ruleHtml).join('')}
              <div style="font-size:0.78rem;color:var(--text-muted);margin-top:10px;margin-bottom:6px;font-weight:600;">REGLAS DE VEBER (2002)</div>
              ${d.veber_rules.map(ruleHtml).join('')}
              <div class="verdict-box ${d.drug_like ? 'pass' : 'fail'}">${d.verdict}</div>
            </div>`;
    } catch (e) {
        el.innerHTML = errorBox(e);
    }
}

// ─── Tab: ADMET ──────────────────────────────────────────────────────────────
async function loadADMET(mol) {
    const el = document.getElementById('admetContent');
    el.innerHTML = spinner();
    try {
        const res  = await postAPI('/analyze/admet', { smiles: mol.smiles });
        const d    = await res.json();

        const COLORS = ['#4facfe','#b794f4','#f6ad55','#68d391','#fc8181'];
        const labels = d.admet.map(a => a.axis);
        const scores = d.admet.map(a => a.score);

        el.innerHTML = `
            <div style="text-align:center;margin-bottom:12px;">
              <span style="font-size:1.8rem;font-weight:800;background:linear-gradient(135deg,#4facfe,#b794f4);-webkit-background-clip:text;-webkit-text-fill-color:transparent;">${d.overall_score}</span>
              <span style="font-size:0.8rem;color:var(--text-muted);display:block;margin-top:2px;">${d.interpretation}</span>
            </div>
            <div class="canvas-wrap"><canvas id="admetChart"></canvas></div>
            <div class="admet-legend">
              ${d.admet.map((a,i) => `
                <div class="admet-item">
                  <div class="admet-item-name" style="color:${COLORS[i]}">${a.axis}</div>
                  <div class="admet-bar-track">
                    <div class="admet-bar-fill" style="width:${a.score}%;background:${COLORS[i]}"></div>
                  </div>
                  <div class="admet-item-desc">${a.description}</div>
                </div>`).join('')}
            </div>`;

        // Destruir gráfico anterior
        if (radarChart) { radarChart.destroy(); radarChart = null; }

        const ctx = document.getElementById('admetChart').getContext('2d');
        radarChart = new Chart(ctx, {
            type: 'radar',
            data: {
                labels,
                datasets: [{
                    label: mol.name,
                    data: scores,
                    backgroundColor: 'rgba(79,172,254,0.15)',
                    borderColor: '#4facfe',
                    pointBackgroundColor: COLORS,
                    pointRadius: 5,
                    borderWidth: 2,
                }]
            },
            options: {
                responsive: true,
                scales: {
                    r: {
                        min: 0, max: 100,
                        ticks: { stepSize: 25, color: '#64748b', font: { size: 10 } },
                        grid:  { color: 'rgba(255,255,255,0.07)' },
                        angleLines: { color: 'rgba(255,255,255,0.07)' },
                        pointLabels: { color: '#94a3b8', font: { size: 12 } },
                    }
                },
                plugins: { legend: { display: false } },
            }
        });
    } catch (e) {
        el.innerHTML = errorBox(e);
    }
}

// ─── Tab: Similitud ──────────────────────────────────────────────────────────
async function loadSimilarity(mol) {
    const el = document.getElementById('similarityContent');
    el.innerHTML = spinner();
    try {
        const res  = await postAPI('/analyze/similarity', { smiles: mol.smiles, top_n: 12 });
        const d    = await res.json();

        if (!d.similar_molecules.length) {
            el.innerHTML = '<div style="color:var(--text-muted);padding:20px;">Sin resultados.</div>';
            return;
        }

        el.innerHTML = `
            <div style="font-size:0.75rem;color:var(--text-muted);margin-bottom:10px;">
              Top ${d.similar_molecules.length} moléculas similares por fingerprint Morgan (radio 2, 2048 bits) + similitud Tanimoto
            </div>
            ${d.similar_molecules.map(m => `
                <div class="sim-card" onclick="selectMolecule(${m.id})">
                  <div class="sim-score" style="color:${simColor(m.similarity)}">${(m.similarity*100).toFixed(1)}%</div>
                  <div class="sim-bar-wrap">
                    <div class="sim-name">${m.name}</div>
                    <div class="sim-track">
                      <div class="sim-fill" style="width:${m.similarity*100}%;background:${simGrad(m.similarity)}"></div>
                    </div>
                  </div>
                  <div class="sim-sol-badge" style="background:${hexToRgba(m.solubility_color,0.2)};color:${m.solubility_color};border:1px solid ${hexToRgba(m.solubility_color,0.4)};">${m.solubility_class}</div>
                </div>`).join('')}`;
    } catch (e) {
        el.innerHTML = errorBox(e);
    }
}

function simColor(s) { return s >= 0.7 ? '#68d391' : s >= 0.4 ? '#f6ad55' : '#94a3b8'; }
function simGrad(s)  { return s >= 0.7 ? 'linear-gradient(90deg,#43e97b,#38f9d7)' : 'linear-gradient(90deg,#4facfe,#00f2fe)'; }

// ─── Tab: Descriptores ───────────────────────────────────────────────────────
async function loadDescriptors(mol) {
    const el = document.getElementById('descriptorsContent');
    el.innerHTML = spinner();
    try {
        const res  = await postAPI('/analyze/descriptors', { smiles: mol.smiles });
        const d    = await res.json();
        allDescriptors = d.descriptors;
        renderDescTable(allDescriptors);
    } catch (e) {
        el.innerHTML = errorBox(e);
    }
}

function renderDescTable(descs) {
    const el = document.getElementById('descriptorsContent');
    el.innerHTML = `
        <div style="font-size:0.75rem;color:var(--text-muted);margin-bottom:8px;">${descs.length} descriptores</div>
        <div class="desc-table-wrap">
          <table>
            <thead><tr><th>Descriptor</th><th>Valor</th></tr></thead>
            <tbody>
              ${descs.map(d => `
                <tr>
                  <td class="td-name">${d.name}</td>
                  <td class="td-val">${d.value !== null ? d.value : '<span style="color:var(--text-muted)">N/A</span>'}</td>
                </tr>`).join('')}
            </tbody>
          </table>
        </div>`;
}

function setupDescSearch() {
    document.getElementById('descSearch').addEventListener('input', e => {
        const q = e.target.value.toLowerCase();
        const filtered = allDescriptors.filter(d => d.name.toLowerCase().includes(q));
        renderDescTable(filtered);
    });
}

// ─── Eventos principales ──────────────────────────────────────────────────────

document.getElementById('imgDownloadBtn')?.addEventListener('click', downloadStructureImage);

document.getElementById('btnDeleteMol')?.addEventListener('click', async () => {
    if (!currentMolecule) return;
    
    const key = prompt("Ingrese la clave para eliminar la molécula:");
    if (key !== "Holamundo") {
        showToast("Clave incorrecta", true);
        return;
    }

    const confirmDelete = confirm(`¿Estás seguro de que deseas eliminar permanentemente '${currentMolecule.name}' de la biblioteca?`);
    if (!confirmDelete) return;

    try {
        const res = await fetch(`${API}/molecule/delete`, {
            method: 'DELETE',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ smiles: currentMolecule.smiles })
        });
        
        const data = await res.json();
        
        if (data.success) {
            showToast(`🗑️ ${currentMolecule.name} eliminada.`);
            document.getElementById('analysisPanel').style.display = 'none';
            document.getElementById('welcomeScreen').style.display = 'flex';
            currentMolecule = null;
            await loadMoleculeLibrary();
        } else {
            showToast(data.error || "Error al eliminar la molécula", true);
        }
    } catch (e) {
        showToast("Error de conexión al eliminar", true);
    }
});

// ─── Tab: Subestructura ──────────────────────────────────────────────────────
function setupSubstructure() {
    document.getElementById('smartsSearchBtn').addEventListener('click', runSubstructure);
    document.getElementById('smartsInput').addEventListener('keydown', e => {
        if (e.key === 'Enter') runSubstructure();
    });
    document.querySelectorAll('.preset-btn').forEach(btn => {
        btn.addEventListener('click', () => {
            document.getElementById('smartsInput').value = btn.dataset.smarts;
            runSubstructure();
        });
    });
}

async function runSubstructure() {
    const smarts = document.getElementById('smartsInput').value.trim();
    if (!smarts) return;
    const el = document.getElementById('substructureContent');
    el.innerHTML = spinner();
    try {
        const res  = await postAPI('/analyze/substructure', {
            smiles: currentMolecule?.smiles || '',
            smarts
        });
        const d    = await res.json();
        if (!d.matches.length) {
            el.innerHTML = `<div style="color:var(--text-muted);padding:12px;font-size:0.83rem;">Sin coincidencias para <code>${smarts}</code></div>`;
            return;
        }
        el.innerHTML = `
            <div style="font-size:0.78rem;color:var(--text-muted);margin:10px 0 6px;">
              ${d.count} coincidencia${d.count !== 1 ? 's' : ''} para <code style="color:var(--accent-teal)">${smarts}</code>
            </div>
            ${d.matches.map(m => `
                <div class="substruct-match" onclick="selectMolecule(${m.id})">
                  <div class="sol-dot" style="background:${m.solubility_color}"></div>
                  <div style="flex:1;font-weight:500;">${m.name}</div>
                  <div style="font-size:0.65rem;padding:2px 8px;border-radius:10px;background:${hexToRgba(m.solubility_color,0.2)};color:${m.solubility_color}">${m.solubility_class}</div>
                </div>`).join('')}`;
    } catch (e) {
        el.innerHTML = errorBox(e);
    }
}

// ─── Helpers ─────────────────────────────────────────────────────────────────
async function postAPI(endpoint, body) {
    return fetch(`${API}${endpoint}`, {
        method:  'POST',
        headers: { 'Content-Type': 'application/json' },
        body:    JSON.stringify(body),
    });
}

function spinner() {
    return '<div class="spinner-wrap"><div class="spinner"></div><span style="color:var(--text-muted);font-size:0.8rem;">Calculando…</span></div>';
}

function errorBox(err) {
    return `<div class="error-box">⚠️ ${err?.message || 'Error desconocido'}</div>`;
}

function hexToRgba(hex, alpha) {
    if (!hex || !hex.startsWith('#')) return `rgba(99,179,237,${alpha})`;
    const r = parseInt(hex.slice(1,3),16);
    const g = parseInt(hex.slice(3,5),16);
    const b = parseInt(hex.slice(5,7),16);
    return `rgba(${r},${g},${b},${alpha})`;
}

// ─── Molecule Builder ─────────────────────────────────────────────────────────

let _validatedSmiles = null; // SMILES validado listo para guardar

function setupMoleculeBuilder() {
    const fab     = document.getElementById('fabNewMol');
    const overlay = document.getElementById('molBuilderOverlay');
    const closeBtn = document.getElementById('modalClose');
    const validateBtn = document.getElementById('validateBtn');
    const saveBtn = document.getElementById('saveBtn');

    fab.addEventListener('click', openBuilderModal);
    closeBtn.addEventListener('click', closeBuilderModal);
    overlay.addEventListener('click', e => { if (e.target === overlay) closeBuilderModal(); });
    document.addEventListener('keydown', e => { if (e.key === 'Escape') closeBuilderModal(); });

    validateBtn.addEventListener('click', validateMolecule);
    saveBtn.addEventListener('click', saveMolecule);

    // Allow Enter in SMILES field to trigger validate
    document.getElementById('molSmilesInput').addEventListener('keydown', e => {
        if (e.key === 'Enter') validateMolecule();
    });
}

function openBuilderModal(prefilledSmiles = '') {
    _validatedSmiles = null;
    document.getElementById('molNameInput').value  = '';
    document.getElementById('molSmilesInput').value = typeof prefilledSmiles === 'string' ? prefilledSmiles : '';
    document.getElementById('molLogsInput').value  = '';
    document.getElementById('validationResult').style.display = 'none';
    document.getElementById('validationResult').innerHTML = '';
    document.getElementById('saveBtn').style.display = 'none';
    document.getElementById('molBuilderOverlay').classList.add('open');
    setTimeout(() => document.getElementById('molNameInput').focus(), 150);
}

function closeBuilderModal() {
    document.getElementById('molBuilderOverlay').classList.remove('open');
}

async function validateMolecule() {
    const name   = document.getElementById('molNameInput').value.trim();
    const smiles = document.getElementById('molSmilesInput').value.trim();
    const resultEl = document.getElementById('validationResult');
    const saveBtn  = document.getElementById('saveBtn');

    if (!smiles) {
        showToast('Ingresa un SMILES para validar.', true);
        document.getElementById('molSmilesInput').focus();
        return;
    }

    resultEl.style.display = 'block';
    resultEl.innerHTML = spinner();
    saveBtn.style.display = 'none';
    _validatedSmiles = null;

    try {
        const res  = await postAPI('/molecule/validate', { smiles, name: name || 'Sin nombre' });
        const data = await res.json();

        if (!data.valid) {
            resultEl.innerHTML = `
                <div class="val-result-invalid">
                    <div class="val-invalid-icon">⚠️</div>
                    <div class="val-invalid-title">Molécula no válida</div>
                    <div class="val-invalid-msg">${data.error || 'RDKit no puede interpretar este SMILES.'}<br><br>
                    Verifica que la cadena SMILES sea correcta. Puedes copiarla desde 
                    <a href="https://pubchem.ncbi.nlm.nih.gov" target="_blank" style="color:var(--accent-blue)">PubChem</a>.</div>
                </div>`;
            return;
        }

        // Molécula válida
        _validatedSmiles = data.smiles; // SMILES canónico de RDKit
        const lip = data.lipinski;
        const sol = data.solubility;
        const p   = data.properties;

        resultEl.innerHTML = `
            <div class="val-result-valid">
                <div class="val-img-wrap">
                    <img src="data:image/png;base64,${data.image_base64}" alt="Estructura 2D">
                </div>
                <div class="val-info">
                    <div style="font-size:0.82rem;font-weight:600;color:var(--text-secondary);margin-bottom:6px;">
                        ✅ Molécula válida — <span style="color:var(--accent-teal);font-family:'Courier New',monospace;font-size:0.78rem;">${data.formula}</span>
                    </div>
                    <div class="val-verdict ${lip.drug_like ? 'ok' : 'warn'}">${lip.verdict}</div>
                    ${sol.logS_esol !== null ? `
                    <div style="font-size:0.75rem;color:var(--text-muted);margin-bottom:8px;">
                        Solubilidad ESOL predicha: 
                        <span style="color:${sol.solubility_color};font-weight:600;">
                            logS = ${sol.logS_esol} (${sol.solubility_class})
                        </span>
                    </div>` : ''}
                    <div class="val-props">
                        <div class="val-prop"><div class="val-prop-label">Peso Mol.</div><div class="val-prop-value">${p.MW}</div></div>
                        <div class="val-prop"><div class="val-prop-label">LogP</div><div class="val-prop-value">${p.LogP}</div></div>
                        <div class="val-prop"><div class="val-prop-label">TPSA (Ų)</div><div class="val-prop-value">${p.TPSA}</div></div>
                        <div class="val-prop"><div class="val-prop-label">HBD</div><div class="val-prop-value">${p.HBD}</div></div>
                        <div class="val-prop"><div class="val-prop-label">HBA</div><div class="val-prop-value">${p.HBA}</div></div>
                        <div class="val-prop"><div class="val-prop-label">Rot. Bonds</div><div class="val-prop-value">${p.RotBonds}</div></div>
                    </div>
                </div>
            </div>`;

        saveBtn.style.display = 'flex';
    } catch (e) {
        resultEl.innerHTML = `<div class="error-box">⚠️ Error al comunicarse con el servidor Python: ${e.message}</div>`;
    }
}

async function saveMolecule() {
    if (!_validatedSmiles) {
        showToast('Valida la molécula primero.', true);
        return;
    }

    const key = prompt("Ingrese la clave para agregar la molécula:");
    if (key !== "Holamundo") {
        showToast("Clave incorrecta", true);
        return;
    }

    const name   = document.getElementById('molNameInput').value.trim();
    const logsRaw = document.getElementById('molLogsInput').value.trim();
    const saveBtn = document.getElementById('saveBtn');

    if (!name) {
        showToast('Escribe un nombre para la molécula.', true);
        document.getElementById('molNameInput').focus();
        return;
    }

    saveBtn.disabled = true;
    saveBtn.textContent = '⏳ Guardando…';

    try {
        const body = { smiles: _validatedSmiles, name };
        if (logsRaw !== '') body.logs_measured = parseFloat(logsRaw);

        const res  = await postAPI('/molecule/add', body);
        const data = await res.json();

        if (data.error) {
            showToast(`❌ ${data.error}`, true);
        } else {
            // Agregar localmente a la lista y refrescar
            if (data.molecule) {
                allMolecules.push(data.molecule);
                renderMolList(allMolecules);
                updateGlobalStats(allMolecules);
                document.getElementById('dbStatus').textContent = `✅ ${allMolecules.length} compuestos cargados`;
            }
            closeBuilderModal();
            showToast(`✅ "${name}" agregada a la biblioteca`);
        }
    } catch (e) {
        showToast('❌ Error de conexión con Python', true);
    } finally {
        saveBtn.disabled = false;
        saveBtn.innerHTML = '💾 Guardar en Biblioteca';
    }
}

// ─── Toast ───────────────────────────────────────────────────────────────────
let _toastTimer = null;
function showToast(msg, isError = false) {
    const toast = document.getElementById('toast');
    toast.textContent = msg;
    toast.classList.toggle('error', isError);
    toast.classList.add('show');
    clearTimeout(_toastTimer);
    _toastTimer = setTimeout(() => toast.classList.remove('show'), 3500);
}

// ─── Reaction Simulator ───────────────────────────────────────────────────────

// Metadata loaded from /api/reactions
let _rxnMeta = {};

async function loadReactionMeta() {
    try {
        const res = await fetch(`${API}/reactions`);
        const list = await res.json();
        list.forEach(r => { _rxnMeta[r.key] = r; });
    } catch (e) {
        // Silently fail — metadata is optional, we fallback gracefully
    }
}

function setupReactionSimulator() {
    const fab      = document.getElementById('fabReaction');
    const overlay  = document.getElementById('reactionOverlay');
    const closeBtn = document.getElementById('rxnModalClose');
    const select   = document.getElementById('rxnTypeSelect');
    const simBtn   = document.getElementById('rxnSimulateBtn');
    const useCurrentBtn = document.getElementById('rxnUseCurrentBtn');

    fab.addEventListener('click', openReactionModal);
    closeBtn.addEventListener('click', closeReactionModal);
    overlay.addEventListener('click', e => { if (e.target === overlay) closeReactionModal(); });
    document.addEventListener('keydown', e => { if (e.key === 'Escape') closeReactionModal(); });

    select.addEventListener('change', updateReactionUI);
    simBtn.addEventListener('click', runReaction);

    useCurrentBtn.addEventListener('click', () => {
        if (currentMolecule) {
            document.getElementById('rxnR1').value = currentMolecule.smiles;
            showToast(`📌 "${currentMolecule.name}" cargado como Reactivo 1`);
        } else {
            showToast('Selecciona primero una molécula del panel lateral.', true);
        }
    });

    // Handle clicks on "Save Product" buttons inside the result area
    document.getElementById('rxnResult').addEventListener('click', e => {
        const saveBtn = e.target.closest('.rxn-save-product');
        if (saveBtn) {
            const smiles = saveBtn.dataset.smiles;
            closeReactionModal();
            openBuilderModal(smiles);
        }
    });
}

function openReactionModal() {
    document.getElementById('rxnR1').value = '';
    document.getElementById('rxnR2').value = '';
    document.getElementById('rxnResult').innerHTML = '';
    // Pre-fill R1 with current molecule if available
    if (currentMolecule) {
        document.getElementById('rxnR1').value = currentMolecule.smiles;
    }
    updateReactionUI();
    document.getElementById('reactionOverlay').classList.add('open');
    setTimeout(() => document.getElementById('rxnR1').focus(), 150);
}

function closeReactionModal() {
    document.getElementById('reactionOverlay').classList.remove('open');
}

function updateReactionUI() {
    const key = document.getElementById('rxnTypeSelect').value;
    const meta = _rxnMeta[key] || {};

    // Update explanation
    const explEl = document.getElementById('rxnExplanation');
    explEl.textContent = meta.explanation || '';

    // Update hints
    document.getElementById('rxnR1Hint').textContent = meta.r1_hint ? `📋 ${meta.r1_hint}` : '';
    document.getElementById('rxnR2Hint').textContent = meta.r2_hint ? `📋 ${meta.r2_hint}` : '';

    // Show/hide R2 and plus
    const needsTwo = meta.needs_two !== false; // default true if missing
    const r2Field  = document.getElementById('rxnR2Field');
    const plus     = document.getElementById('rxnPlus');
    r2Field.style.display = needsTwo ? '' : 'none';
    plus.style.display    = needsTwo ? '' : 'none';
}

async function runReaction() {
    const rxnKey  = document.getElementById('rxnTypeSelect').value;
    const smiles1 = document.getElementById('rxnR1').value.trim();
    const smiles2 = document.getElementById('rxnR2').value.trim();
    const resultEl = document.getElementById('rxnResult');

    if (!smiles1) {
        showToast('Ingresa el SMILES del Reactivo 1.', true);
        document.getElementById('rxnR1').focus();
        return;
    }

    const meta = _rxnMeta[rxnKey] || {};
    if (meta.needs_two !== false && !smiles2) {
        showToast('Esta reacción requiere Reactivo 2.', true);
        document.getElementById('rxnR2').focus();
        return;
    }

    resultEl.innerHTML = spinner();

    try {
        const res  = await postAPI('/analyze/reaction', {
            reaction_type: rxnKey,
            reactant1: smiles1,
            reactant2: smiles2,
        });
        const data = await res.json();

        if (data.error) {
            resultEl.innerHTML = `
                <div class="error-box">
                    ⚠️ ${data.error}
                    ${data.hint ? `<div class="rxn-error-hint">${data.hint}</div>` : ''}
                </div>`;
            return;
        }

        renderReactionResult(data);
    } catch (e) {
        resultEl.innerHTML = errorBox(e);
    }
}

function renderReactionResult(data) {
    const resultEl = document.getElementById('rxnResult');

    const reactantCard = (r, label) => r ? `
        <div class="rxn-mol-card reactant">
            <div class="rxn-mol-label">${label}</div>
            <div class="rxn-mol-formula">${r.formula || ''}</div>
            ${r.image_base64
                ? `<img class="rxn-mol-img" src="data:image/png;base64,${r.image_base64}" alt="${label}">`
                : '<div class="rxn-mol-no-img">Sin imagen</div>'}
            <div class="rxn-mol-smiles">${r.smiles || ''}</div>
        </div>` : '';

    const r2Html = data.reactant2
        ? `<div class="rxn-plus-sign">+</div>${reactantCard(data.reactant2, 'Reactivo 2')}`
        : '';

    const productCard = p => `
        <div class="rxn-mol-card product">
            <div class="rxn-mol-label">🟢 Producto</div>
            <div class="rxn-mol-formula">${p.formula}</div>
            ${p.image_base64
                ? `<img class="rxn-mol-img" src="data:image/png;base64,${p.image_base64}" alt="Producto">`
                : '<div class="rxn-mol-no-img">Sin imagen</div>'}
            <div class="rxn-mol-smiles">${p.smiles}</div>
            <div class="rxn-prod-verdict ${p.lipinski.drug_like ? 'ok' : 'warn'}">${p.lipinski.verdict}</div>
            ${p.solubility.logS_esol !== null
                ? `<div class="rxn-prod-sol" style="color:${p.solubility.color}">logS = ${p.solubility.logS_esol} · ${p.solubility.class}</div>`
                : ''}
            <div class="rxn-prod-props">
                <span>MW ${p.properties.MW}</span>
                <span>LogP ${p.properties.LogP}</span>
                <span>TPSA ${p.properties.TPSA}</span>
                <span>HBD ${p.properties.HBD}</span>
                <span>HBA ${p.properties.HBA}</span>
            </div>
            <button class="btn btn-ghost rxn-save-product" data-smiles="${p.smiles}" style="margin-top:10px;font-size:0.75rem;padding:6px 12px;width:100%;justify-content:center;">
                💾 Guardar como Molécula
            </button>
        </div>`;

    resultEl.innerHTML = `
        <div class="rxn-result-header">
            <span class="rxn-result-title">⚗️ ${data.reaction_name}</span>
            <span class="rxn-result-count">${data.products.length} producto${data.products.length !== 1 ? 's' : ''} generado${data.products.length !== 1 ? 's' : ''}</span>
        </div>
        <div class="rxn-explanation-box">${data.explanation}</div>

        ${data.products.map(p => `
        <div class="rxn-scheme">
            ${reactantCard(data.reactant1, 'Reactivo 1')}
            ${r2Html}
            <div class="rxn-arrow">
                <div class="rxn-arrow-line"></div>
                <div class="rxn-arrow-head">▶</div>
            </div>
            ${productCard(p)}
        </div>`).join('')}
    `;
}

// ─── SHUTDOWN APP ────────────────────────────────────────────────────────────
document.getElementById('btnExitApp')?.addEventListener('click', async () => {
    if (!confirm("¿Deseas apagar los servidores y salir de la aplicación?")) return;
    try {
        await fetch('/shutdown', { method: 'POST' });
        document.body.innerHTML = '<div style="display:flex;height:100vh;align-items:center;justify-content:center;background:#0f172a;color:white;font-family:Inter,sans-serif;flex-direction:column;"><h1>Servidores apagados exitosamente.</h1><p style="color:#94a3b8;margin-top:10px;">Ya puedes cerrar esta pestaña.</p></div>';
    } catch(e) {
        showToast("Se envió la señal de apagado", true);
        setTimeout(() => {
            document.body.innerHTML = '<div style="display:flex;height:100vh;align-items:center;justify-content:center;background:#0f172a;color:white;font-family:Inter,sans-serif;flex-direction:column;"><h1>Servidores apagados exitosamente.</h1><p style="color:#94a3b8;margin-top:10px;">Ya puedes cerrar esta pestaña.</p></div>';
        }, 1000);
    }
});
