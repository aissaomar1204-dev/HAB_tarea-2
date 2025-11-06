# Network Propagation - Aissa Omar

## Descripción

Este proyecto implementa algoritmos de **propagación en redes** para la priorización de genes asociados a enfermedades. Se han implementado dos algoritmos principales:

1. **GUILD** (Genes Underlying Inheritance Linked Disorders) - usando Random Walk with Restart
2. **DIAMOnD** (DIseAse MOdule Detection) - expansión iterativa de módulos

## Genes Semilla

Los genes utilizados como semillas son:
- **ENO1** (Enolase 1)
- **PGK1** (Phosphoglycerate Kinase 1)  
- **HK2** (Hexokinase 2)

Estos genes están relacionados con la vía metabólica de la glucólisis.

## Estructura del Proyecto

```
network_propagation/
├── data/
│   ├── network_guild.txt                        # Red formateada para GUILD
│   ├── network_diamond.txt                      # Red formateada para DIAMOnD
│   ├── string_network_filtered_hugo-400.tsv     # Red filtrada de STRING
│   └── genes_seed.txt                           # Genes semilla
├── scripts/
│   ├── process_STRING.py                        # Procesamiento de red STRING
│   └── script.py                             # Script principal ⭐
├── results/                                     # Resultados generados
│   ├── guild_results.txt
│   └── diamond_results.txt
├── README.md
└── requirements.txt
```

## Instalación

### Requisitos Previos
- Python 3.7 o superior
- pip

### Instalación de Dependencias

```bash
pip install -r requirements.txt
```

## Uso

### Sintaxis General

```bash
python scripts/script.py --algorithm [guild|diamond|both] \
                            --network <archivo_red> \
                            --seeds <archivo_semillas> \
                            --output <archivo_salida>
```

### Ejemplos de Uso

#### 1. Ejecutar solo GUILD

```bash
python scripts/script.py \
    --algorithm guild \
    --network data/network_guild.txt \
    --seeds data/genes_seed.txt \
    --output results/guild_results.txt \
    --alpha 0.85
```

#### 2. Ejecutar solo DIAMOnD

```bash
python scripts/script.py \
    --algorithm diamond \
    --network data/network_diamond.txt \
    --seeds data/genes_seed.txt \
    --output results/diamond_results.txt \
    --max-genes 200
```

#### 3. Ejecutar ambos algoritmos

```bash
python scripts/script.py \
    --algorithm both \
    --network data/network_guild.txt \
    --seeds data/genes_seed.txt \
    --output results/results.txt
```

### Parámetros Disponibles

| Parámetro | Descripción | Default | Requerido |
|-----------|-------------|---------|-----------|
| `--algorithm` | Algoritmo a ejecutar: `guild`, `diamond`, o `both` | `both` | No |
| `--network` | Archivo de red de interacciones | - | **Sí** |
| `--seeds` | Archivo con genes semilla (uno por línea) | - | **Sí** |
| `--output` | Archivo de salida para resultados | - | **Sí** |
| `--alpha` | Parámetro de restart para GUILD (0-1) | 0.85 | No |
| `--max-genes` | Número máximo de genes para DIAMOnD | 200 | No |
| `--network-type` | Tipo de formato: `auto`, `guild`, `diamond`, `string` | `auto` | No |

## Formatos de Archivos

### Formato de Red

El archivo de red debe tener tres columnas separadas por espacios o tabuladores:

```
node1   node2   weight
ENO1    PGK1    0.95
PGK1    HK2     0.87
...
```

### Formato de Semillas

El archivo de semillas debe contener un gen por línea:

```
ENO1
PGK1
HK2
```

### Formato de Salida

#### GUILD
```
Rank    Gene    GUILD_Score
1       ENO1    0.15234
2       PGK1    0.14876
3       GAPDH   0.08932
...
```

#### DIAMOnD
```
Rank    Gene    DIAMOnD_Iteration    Connectivity
1       ENO1    0                    inf
2       PGK1    0                    inf
3       HK2     0                    inf
4       GAPDH   1                    2.85
...
```

## Algoritmos Implementados

### 1. GUILD - Random Walk with Restart (RWR)

**Descripción**: Simula un caminante aleatorio en la red que puede moverse a nodos vecinos o reiniciar en los nodos semilla.

**Ecuación matemática**:
```
p(t+1) = α × W × p(t) + (1-α) × p(0)
```

Donde:
- `p(t)`: Vector de probabilidades en el tiempo t
- `W`: Matriz de transición normalizada por filas
- `p(0)`: Vector inicial (distribución uniforme en semillas)
- `α`: Parámetro de restart (típicamente 0.7-0.9)

**Ventajas**:
- Captura la estructura global de la red
- Prioriza genes cercanos a las semillas
- Robusto ante ruido en la red

**Referencias**:
- Köhler et al. (2008) "Walking the Interactome for Prioritization of Candidate Disease Genes"
- Guney et al. (2016) "Network-based in silico drug efficacy screening"

### 2. DIAMOnD - DIseAse MOdule Detection

**Descripción**: Algoritmo iterativo que expande un módulo de enfermedad añadiendo genes con alta conectividad.

**Algoritmo**:
1. Iniciar con genes semilla S
2. Para cada iteración:
   - Calcular conectividad de cada gen candidato al módulo
   - Añadir el gen con mayor conectividad
3. Repetir hasta alcanzar el tamaño deseado

**Ventajas**:
- Identifica módulos de enfermedad coherentes
- Preserva la estructura local de la red
- Interpretación biológica clara

**Referencias**:
- Ghiassian et al. (2015) "A DIseAse MOdule Detection (DIAMOnD) Algorithm Derived from a Systematic Analysis of Connectivity Patterns of Disease Proteins in the Human Interactome"

## Interpretación de Resultados

### Scores GUILD
- Valores más altos indican mayor relevancia funcional
- Los genes semilla típicamente tienen los scores más altos
- Genes con scores >0.01 son candidatos relevantes

### DIAMOnD Iterations
- Iteración 0: Genes semilla
- Iteraciones bajas (1-50): Alta confianza
- Iteraciones medias (51-100): Confianza moderada
- Iteraciones altas (>100): Exploración de la red

### Connectivity
- Suma de pesos de aristas entre el gen y el módulo
- Valores más altos indican mayor integración con el módulo
- Útil para filtrar candidatos por relevancia

## Librerías Utilizadas

- **NetworkX** (v2.x): Análisis y manipulación de redes
  - Creación de grafos
  - Cálculo de matriz de adyacencia
  - Operaciones de red

- **Pandas** (v1.x): Manipulación de datos tabulares
  - Lectura/escritura de archivos
  - Procesamiento de resultados
  - Formateo de salida

- **NumPy** (v1.x): Cálculo numérico
  - Operaciones matriciales
  - Álgebra lineal (para RWR)
  - Vectorización de operaciones

## Conceptos Clave

### Network Propagation
Técnica que explota la estructura de redes biológicas para propagar información desde genes conocidos (semillas) hacia genes candidatos relacionados.

### Random Walk
Proceso estocástico que modela el movimiento aleatorio en un grafo, útil para medir proximidad funcional entre nodos.

### Disease Module
Conjunto de genes funcionalmente relacionados que contribuyen a una enfermedad. Típicamente forman subgrafos densos en la red de interacciones.

## ⚠️ Notas Importantes

1. **Calidad de la Red**: Los resultados dependen críticamente de la calidad y completitud de la red de interacciones.

2. **Selección de Parámetros**:
   - `alpha` en GUILD: Valores típicos 0.7-0.9. Valores altos favorecen exploración global.
   - `max_genes` en DIAMOnD: Depende del tamaño esperado del módulo (típicamente 100-500).

3. **Tiempo de Ejecución**: 
   - GUILD es más rápido (convergencia iterativa)
   - DIAMOnD puede ser más lento en redes grandes

4. **Semillas No Encontradas**: Si alguna semilla no está en la red, se mostrará una advertencia pero el análisis continuará.

## Troubleshooting

### Error: "No hay genes semilla cargados"
**Solución**: Verificar que el archivo de semillas existe y contiene genes válidos.

### Error: "Ninguna semilla encontrada en la red"
**Solución**: Verificar que los nombres de genes en semillas coinciden con los de la red.

### Warning: "No se alcanzó convergencia"
**Solución**: Aumentar `max_iter` o disminuir `tol` en el código.

## Autor

Aissa Omar El Hammoouti Chachoui

## Referencias

1. Ghiassian, S. D., et al. (2015). "A DIseAse MOdule Detection (DIAMOnD) Algorithm..." PLOS Computational Biology.

2. Guney, E., et al. (2016). "Network-based in silico drug efficacy screening." Nature Communications.

3. Köhler, S., et al. (2008). "Walking the Interactome for Prioritization of Candidate Disease Genes." The American Journal of Human Genetics.

4. STRING Database: https://string-db.org/
