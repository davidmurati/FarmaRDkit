/**
 * Node.js/Express Server — RDKit Pharma Dashboard
 * Puerto: 3000 (Frontend) | Proxy a Python Flask :5000
 */
const express = require('express');
const { createProxyMiddleware } = require('http-proxy-middleware');
const path = require('path');

const app = express();
const PORT = 3000;
const PYTHON_API = 'http://localhost:5000';

// Servir archivos estáticos del frontend
app.use(express.static(path.join(__dirname, 'public')));

// Proxy: reenviar /api/* al servidor Python
app.use('/api', createProxyMiddleware({
    target: PYTHON_API,
    changeOrigin: true,
    pathRewrite: { '^/api': '' },
    on: {
        error: (err, req, res) => {
            console.error('[Proxy Error]', err.message);
            res.status(502).json({
                error: 'No se puede conectar con el backend Python.',
                detail: 'Asegúrate de que el servidor Python esté corriendo en el puerto 5000.',
                hint: 'Ejecuta: python python/server.py'
            });
        }
    }
}));

// Ruta principal — servir index.html
app.get('/', (req, res) => {
    res.sendFile(path.join(__dirname, 'public', 'index.html'));
});

app.listen(PORT, () => {
    console.log(`\n🧪 RDKit Pharma Dashboard`);
    console.log(`─────────────────────────────────────`);
    console.log(`🌐 Frontend: http://localhost:${PORT}`);
    console.log(`🔬 Python API: ${PYTHON_API}`);
    console.log(`─────────────────────────────────────\n`);
});
